#include "mpiwrap.h"
#if defined(__PLUMED_MPIWRAP) && !defined(__PLUMED_HAS_MPI)
// Provider selection from injection, the environment, or native communicator import.
#include <cstdlib>
#include <cstring>
#include <stdexcept>
#include <string>
#include <mutex>
#include <condition_variable>
#include <thread>
#include <vector>
#include "internal.h"
#if defined(__PLUMED_HAS_DLOPEN)
#include "../config/Config.h"
#include <dlfcn.h>
#include <spawn.h>
#include <sys/stat.h>
#include <sys/wait.h>
#include <unistd.h>
#include <cerrno>
#include <memory>
#ifdef __APPLE__
#include <fstream>
#include <mach-o/loader.h>
#include <mach-o/fat.h>
#include <libkern/OSByteOrder.h>
#endif
extern char **environ;
#endif

namespace PLMD {
namespace mpiwrap {
namespace {

#if defined(__PLUMED_HAS_DLOPEN)
// dyld can open an executable compiler launcher as though it were a library.
// Inspect its file type before dlopen, including universal Mach-O containers.
bool isCompilerExecutable(const char *path) {
#ifdef __APPLE__
  std::ifstream input(path, std::ios::binary);
  const auto executable_at = [&input](std::uint64_t offset) {
    mach_header header{};
    input.clear();
    input.seekg(static_cast<std::streamoff>(offset));
    if (!input.read(reinterpret_cast<char *>(&header), sizeof(header))) {
      return false;
    }
    if (header.magic == MH_MAGIC || header.magic == MH_MAGIC_64) {
      return header.filetype == MH_EXECUTE;
    }
    if (header.magic == MH_CIGAM || header.magic == MH_CIGAM_64) {
      return OSSwapInt32(header.filetype) == MH_EXECUTE;
    }
    return false;
  };
  if (executable_at(0)) {
    return true;
  }
  fat_header header{};
  input.clear();
  input.seekg(0);
  if (!input.read(reinterpret_cast<char *>(&header), sizeof(header))) {
    return false;
  }
  const bool swapped = header.magic == FAT_CIGAM || header.magic == FAT_CIGAM_64;
  const bool wide = header.magic == FAT_MAGIC_64 || header.magic == FAT_CIGAM_64;
  if (!swapped && header.magic != FAT_MAGIC && header.magic != FAT_MAGIC_64) {
    return false;
  }
  const std::uint32_t count = swapped ? OSSwapInt32(header.nfat_arch) : header.nfat_arch;
  if (!count || count > 128) {
    return false;
  }
  for (std::uint32_t i = 0; i < count; ++i) {
    input.clear();
    input.seekg(sizeof(header) + i * (wide ? sizeof(fat_arch_64) : sizeof(fat_arch)));
    std::uint64_t offset;
    if (wide) {
      fat_arch_64 arch{};
      if (!input.read(reinterpret_cast<char *>(&arch), sizeof(arch))) {
        return false;
      }
      offset = swapped ? OSSwapInt64(arch.offset) : arch.offset;
    } else {
      fat_arch arch{};
      if (!input.read(reinterpret_cast<char *>(&arch), sizeof(arch))) {
        return false;
      }
      offset = swapped ? OSSwapInt32(arch.offset) : arch.offset;
    }
    if (!executable_at(offset)) {
      return false;
    }
  }
  return true;
#else
  (void) path;
  return false;
#endif
}

struct ProviderBuildDirectory {
  std::string directory;
  std::string library;
  ProviderBuildDirectory() {
    const char *temporary_root = std::getenv("TMPDIR");
    std::string pattern = std::string(temporary_root && temporary_root[0] ?
                                      temporary_root : "/tmp") + "/plumed-mpi-provider.XXXXXX";
    // mkdtemp edits its argument; C++17 guarantees writable string storage.
    if (!mkdtemp(&pattern[0])) {
      throw std::runtime_error("could not create MPI provider build directory");
    }
    directory = pattern;
#ifdef __APPLE__
    library = directory + "/provider.dylib";
#else
    library = directory + "/provider.so";
#endif
  }
  ~ProviderBuildDirectory() {
    unlink(library.c_str());
    rmdir(directory.c_str());
  }
};

void installProviderLibrary(Registry &registry, void *library) {
  using Factory = const plumed_mpi_provider *(*)();
  void *symbol = dlsym(library, PLUMED_MPIW_PROVIDER_FACTORY_SYMBOL);
  Factory factory = nullptr;
  static_assert(sizeof(factory) == sizeof(symbol), "dlopen function pointer representation");
  std::memcpy(&factory, &symbol, sizeof(factory));
  try {
    if (!factory) {
      throw std::runtime_error("selected MPI provider library has no provider factory");
    }
    std::string error;
    if (!registry.install(factory(), error, library)) {
      throw std::runtime_error("selected MPI provider is incompatible: " + error);
    }
  } catch (...) {
    dlclose(library);
    throw;
  }
  // Successful installation transfers this dlopen reference to the registry.
}

bool installAdjacentProvider(Registry &registry) {
  const std::string kernel = config::getLibraryPath();
  if (kernel.empty()) {
    return false;
  }
  // Follow a kernel symlink and resolve relative loader paths before looking
  // beside the actual library. Do not depend on the active Conda environment.
  const std::unique_ptr<char, decltype(&std::free)> resolved(
    realpath(kernel.c_str(), nullptr), &std::free);
  if (!resolved) {
    return false;
  }
  const std::string location(resolved.get());
  const std::string path = location.substr(0, location.find_last_of('/') + 1) +
                           "lib" + config::getPlumedProgramName() + "-mpi-provider." + config::getSoExt();
  struct stat information;
  // lstat also detects a broken provider symlink: an installed but unusable
  // provider must report its error rather than silently compiling another one.
  if (lstat(path.c_str(), &information) != 0) {
    if (errno == ENOENT) {
      return false;
    }
    throw std::runtime_error("could not inspect installed MPI provider " + path +
                             ": " + std::strerror(errno));
  }
  void *library = dlopen(path.c_str(), RTLD_NOW | RTLD_LOCAL);
  if (!library) {
    const char *message = dlerror();
    throw std::runtime_error("could not load installed MPI provider " + path +
                             ": " + (message ? message : "unknown loader error"));
  }
  try {
    installProviderLibrary(registry, library);
  } catch (const std::exception &error) {
    throw std::runtime_error("installed MPI provider " + path + ": " + error.what());
  }
  return true;
}

void compileSelectedProvider(const char *compiler, const std::string &output) {
  std::string script = config::getPlumedRoot() + "/scripts/mkmpi.sh";
  // Spawn an interpreter with literal arguments. Neither compiler paths nor
  // output paths are shell commands, even when they contain spaces or '$'.
  std::vector<char *> arguments{
    const_cast<char *>("/bin/bash"), const_cast<char *>(script.c_str())
  };
  // A null compiler deliberately delegates discovery to mkmpi, including
  // MPICC/MPICXX and the host application's module-provided PATH.
  if (compiler && compiler[0]) {
    arguments.push_back(const_cast<char *>("--compiler"));
    arguments.push_back(const_cast<char *>(compiler));
  }
  arguments.push_back(const_cast<char *>("--out"));
  arguments.push_back(const_cast<char *>(output.c_str()));
  arguments.push_back(nullptr);
  // A kernel loaded directly by a host has not run the plumed shell wrappers.
  // Supply their configured paths to this child only; changing the process
  // environment would race with other callers and unrelated application code.
  const std::string root = "PLUMED_ROOT=" + config::getPlumedRoot();
  const std::string includedir = "PLUMED_INCLUDEDIR=" + config::getPlumedIncludedir();
  const std::string program = "PLUMED_PROGRAM_NAME=" + config::getPlumedProgramName();
  std::vector<char *> child_environment;
  for(char **entry = environ; entry && *entry; ++entry) {
    if(std::strncmp(*entry,"PLUMED_ROOT=",12)!=0 &&
        std::strncmp(*entry,"PLUMED_INCLUDEDIR=",18)!=0 &&
        std::strncmp(*entry,"PLUMED_PROGRAM_NAME=",20)!=0) {
      child_environment.push_back(*entry);
    }
  }
  child_environment.push_back(const_cast<char *>(root.c_str()));
  child_environment.push_back(const_cast<char *>(includedir.c_str()));
  child_environment.push_back(const_cast<char *>(program.c_str()));
  child_environment.push_back(nullptr);
  pid_t child;
  const int error = posix_spawn(&child, "/bin/bash", nullptr, nullptr, arguments.data(),
                                child_environment.data());
  if (error != 0) {
    throw std::runtime_error("could not start MPI provider compiler: " +
                             std::string(std::strerror(error)));
  }
  int status = 0;
  pid_t result;
  do {
    result = waitpid(child, &status, 0);
  } while (result < 0 && errno == EINTR);
  if (result < 0 || !WIFEXITED(status) || WEXITSTATUS(status) != 0) {
    throw std::runtime_error("MPI provider compilation failed; see compiler diagnostics above");
  }
}

void installCompiledProvider(Registry &registry, const char *compiler) {
  ProviderBuildDirectory build;
  compileSelectedProvider(compiler, build.library);
  void *library = dlopen(build.library.c_str(), RTLD_NOW | RTLD_LOCAL);
  if (!library) {
    const char *message = dlerror();
    throw std::runtime_error("could not load compiled MPI provider: " +
                             std::string(message ? message : "unknown loader error"));
  }
  installProviderLibrary(registry, library);
}
#endif

void installEnvironmentProvider(Registry &registry, const char *path) {
  if (!path || path[0] != '/') {
    throw std::runtime_error("PLUMED_MPI_WRAPPER requires an absolute library or compiler path");
  }
#if defined(__PLUMED_HAS_DLOPEN)
  const bool executable = isCompilerExecutable(path);
  void *library = executable ? nullptr : dlopen(path, RTLD_NOW | RTLD_LOCAL);
  if (library) {
    installProviderLibrary(registry, library);
    return;
  }
  const char *message = executable ? "selected file is an executable compiler" : dlerror();
  const std::string load_error = message ? message : "unknown loader error";
  const std::string selected(path);
  const std::string basename = selected.substr(selected.find_last_of('/') + 1);
  const bool library_name =
    (basename.size() >= 3 && basename.compare(basename.size() - 3, 3, ".so") == 0) ||
    basename.find(".so.") != std::string::npos ||
    (basename.size() >= 6 && basename.compare(basename.size() - 6, 6, ".dylib") == 0) ||
    (basename.size() >= 4 && basename.compare(basename.size() - 4, 4, ".dll") == 0);
  struct stat information;
  if (library_name || stat(path, &information) != 0 ||
      !S_ISREG(information.st_mode) || access(path, X_OK) != 0) {
    throw std::runtime_error("could not load selected MPI provider: " + load_error);
  }
  installCompiledProvider(registry, path);
#else
  (void) registry;
  throw std::runtime_error("PLUMED_MPI_WRAPPER requires a kernel built with dlopen support");
#endif
}

} // namespace

bool forceEnvironmentProvider() {
  const char *force = std::getenv("PLUMED_MPI_WRAPPER_FORCE");
  if (!force || !force[0] || std::strcmp(force, "0") == 0) {
    return false;
  }
  if (std::strcmp(force, "1") != 0) {
    throw std::runtime_error("PLUMED_MPI_WRAPPER_FORCE must be 0 or 1");
  }
  const char *path = std::getenv("PLUMED_MPI_WRAPPER");
  if (!path || !path[0]) {
    throw std::runtime_error("PLUMED_MPI_WRAPPER_FORCE=1 requires PLUMED_MPI_WRAPPER");
  }
  static std::mutex selection_mutex;
  static std::string forced_path;
  std::lock_guard<std::mutex> lock(selection_mutex);
  Registry &registry = currentRegistry();
  if (!forced_path.empty()) {
    if (forced_path != path) {
      throw std::runtime_error("cannot change PLUMED_MPI_WRAPPER after forced provider installation");
    }
    return true;
  }
  if (registry.installed()) {
    throw std::runtime_error("PLUMED_MPI_WRAPPER_FORCE must be set before MPI provider installation");
  }
  std::string selected_path(path);
  installEnvironmentProvider(registry, path);
  forced_path.swap(selected_path);
  return true;
}

namespace {
void selectProvider(bool build_default) {
  if (forceEnvironmentProvider()) {
    return;
  }
  static std::mutex selection_mutex;
  static std::condition_variable selection_changed;
  static bool selecting = false;
  static std::thread::id selecting_thread;

  Registry &registry = currentRegistry();
  if (registry.installed()) {
    return;
  }
  const char *selected = std::getenv("PLUMED_MPI_WRAPPER");
  if ((!selected || !selected[0]) && !build_default) {
    return;
  }
  {
    std::unique_lock<std::mutex> lock(selection_mutex);
    // Another thread may install a provider while we wait for the lock.
    // cppcheck-suppress identicalConditionAfterEarlyExit
    if (registry.installed()) {
      return;
    }
    while (selecting) {
      if (selecting_thread == std::this_thread::get_id()) {
        throw std::runtime_error("recursive MPI provider selection");
      }
      selection_changed.wait(lock);
      if (registry.installed()) {
        return;
      }
    }
    selecting = true;
    selecting_thread = std::this_thread::get_id();
  }

  try {
    if (selected && selected[0]) {
      installEnvironmentProvider(registry, selected);
    } else {
#if defined(__PLUMED_HAS_DLOPEN)
      if (!installAdjacentProvider(registry)) {
        installCompiledProvider(registry, nullptr);
      }
#else
      throw std::runtime_error("automatic MPI provider building requires dlopen support");
#endif
    }
  } catch (...) {
    std::lock_guard<std::mutex> lock(selection_mutex);
    selecting = false;
    selecting_thread = std::thread::id();
    selection_changed.notify_all();
    throw;
  }
  {
    std::lock_guard<std::mutex> lock(selection_mutex);
    selecting = false;
    selecting_thread = std::thread::id();
  }
  selection_changed.notify_all();
}

} // namespace

void ensureProviderSelected() {
  selectProvider(false);
}

void ensureProviderForCommunicator() {
  selectProvider(true);
}

} // namespace mpiwrap
} // namespace PLMD

#endif
