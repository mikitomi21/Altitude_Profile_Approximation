import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

data = pd.read_csv("Data/stale.txt", sep=" ")
print(data.head(5))