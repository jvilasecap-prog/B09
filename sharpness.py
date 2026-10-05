import numpy as np
import pandas as pd

files = ["airfoils/ms313.dat.txt", 
         "airfoils/ms317.dat.txt", 
         "airfoils/sc20714.dat.txt", 
         "airfoils/whitcomb.dat.txt"]

for file_path in files:
    with open(file_path, "r") as f:
        airfoil_name = f.readline().strip()

    df = pd.read_csv(
        file_path,
        sep=r"\s+",  # deals with whitespace
        skiprows=1,  # skips the first line (airfoil name)
        names=["x", "y"],  
        header=None,  
    )

    df = df[(df["y"] >= 0) & (df["x"] <= 1)] # upper surface
    df = df.drop_duplicates() # upper surface

    # print(df)

    x = [0.0015, 0.0600]
    y = []

    for xi in x:
        a = xi
        b = xi
        lower_bound_found = False
        upper_bound_found = False
        if(df[(df["x"] == xi)].shape[0] == 1):
            y.append(df[(df["x"] == xi)]["y"])
        else:
            while not lower_bound_found:
                a -= 0.001
                if(df[(df["x"] >= a) & (df["x"] <= xi)].shape[0] >= 1):
                    lower_bound_found = True
            while not upper_bound_found:
                b += 0.001
                if(df[(df["x"] >= a) & (df["x"] <= b)].shape[0] >= 2):
                    upper_bound_found = True
            df_linear = df[(df["x"] >= a) & (df["x"] <= b)]
            # print(df_linear)
            y.append(np.interp(xi, df_linear["x"], df_linear["y"]))
            
            
    dy = (y[1] - y[0]) * 100 # [%]

    print(f"{airfoil_name} \t {dy}")



