# models whose output is not a contiguous tensor
import torch


class Strided(torch.nn.Module):
    def forward(self, x):
        h = torch.cat([torch.sin(x), torch.cos(x)], 1)
        return h[:, ::2]  # [sin(x0), cos(x0)] for two inputs, with stride 2


class Expand(torch.nn.Module):
    def forward(self, x):
        return (x * x).sum(1, keepdim=True).expand(1, 3)  # stride 0


torch.jit.script(Strided()).save("strided.ptc")
torch.jit.script(Expand()).save("expand.ptc")
