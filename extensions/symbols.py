"""Native identity symbol names shared by source and plugin generators."""


def module_symbol(name):
    # Preserve the existing ABI for C-compatible names. A separate prefix plus
    # complete hexadecimal encoding makes hyphenated names collision-free.
    if "-" in name:
        return "plumed_encoded_module_" + name.encode("ascii").hex()
    return "plumed_module_" + name
