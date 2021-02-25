import gpytorch
import numpy as np
import torch

print("module loaded")

def initialize(mol_structures):
    print("Received as input :\n", mol_structures)
    print(torch.version())
    return 1