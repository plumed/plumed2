/* A fake CUDA driver for GpuDevice::callerCudaDevice(): built as libcuda.so.1 by the config */
static int device=-1;
static int context;
void fakeSetDevice(int d) {
  device=d;
}
int cuCtxGetCurrent(void** ctx) {
  *ctx = device>=0 ? (void*)&context : (void*)0;
  return 0;
}
int cuCtxGetDevice(int* d) {
  *d=device;
  return 0;
}
