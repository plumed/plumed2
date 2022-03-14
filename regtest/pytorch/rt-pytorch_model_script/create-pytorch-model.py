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
    if x > 0:
        # CV definition
        y = torch.sin(x)
    else:
        y = torch.tan(x)

    return y

input_size = 1

# Compile via scripting
scripted_cv = torch.jit.script( my_torch_cv )

filename='torch_model.pt'
scripted_cv.save(filename)
 
