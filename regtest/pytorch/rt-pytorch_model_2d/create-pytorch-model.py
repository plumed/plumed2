import torch
print(torch.__version__)

def my_torch_cv(x):
    '''
    Here goes the definition of the CV.

    Inputs:
        x (torch.tensor): input, either scalar or 1-D array
    Return:
        y (torch.tensor): collective variable (scalar)
    '''
    # CV definition
    #y1 = torch.sin(x)
    #y2 = torch.cos(x)
    #y = torch.cat((y1,y2))
    y = torch.sin(x)
    return y

input_size = 2

# -- DEFINE INPUT -- 
#random 
#x = torch.rand(input_size, dtype=torch.float32, requires_grad=True).unsqueeze(0)
#or by choosing the value(s) of the array
x = torch.tensor([0.,1.57], dtype=torch.float32, requires_grad=True)

# -- CALCULATE CV -- 
y = my_torch_cv(x)

# -- CALCULATE DERIVATIVES -- 
for yy in y:
    dy = torch.autograd.grad(yy, x, create_graph=True)
    # -- PRINT -- 
    print('CV TEST')
    print('n_input\t: {}'.format(input_size))
    print('x\t: {}'.format(x))
    print('cv\t: {}'.format(yy))
    print('der\t: {}'.format(dy))

# Compile via tracing
traced_cv   = torch.jit.trace ( my_torch_cv, example_inputs=x )
filename='torch_model.pt'
traced_cv.save(filename)

