from pathlib import Path
from decimal import Decimal, getcontext
import shutil
import numpy as np
import pandas as pd
import subprocess
import os
#this script zip folders by chunks

def format_using_decimal(value, precision=4):
    # Set the precision higher to ensure correct conversion
    getcontext().prec = precision + 2
    # Convert the float to a Decimal with exact precision
    decimal_value = Decimal(str(value))
    # Normalize to remove trailing zeros
    formatted_value = decimal_value.quantize(Decimal(1)) if decimal_value == decimal_value.to_integral() else decimal_value.normalize()
    return str(formatted_value)
N=5 #unit cell number
csv_T_root=f"./dataAll/N{N}/csvOutAll/"

number=int((10-0.1)/0.01)
TVals=[0.1+0.01*n for n in range(0,number+1)]
lastFileNum=20
TStrAll=[]
chunk_size = 50
chunks = [TVals[i:i + chunk_size] for i in range(0, len(TVals), chunk_size)]

for chunk_index, chunk in enumerate(chunks):
    # Convert each T value in the chunk to its formatted string
    chunk_TStrs = [format_using_decimal(T) for T in chunk]
    folder_names = [f"T{TStr}" for TStr in chunk_TStrs]
    # Check for existing folders
    existing_folders = []
    for folder in folder_names:
        folder_path = os.path.join(csv_T_root, folder)
        if os.path.isdir(folder_path):
            existing_folders.append(folder)
        else:
            print(f"Warning: Folder {folder} does not exist and will be skipped.")

    if not existing_folders:
        print(f"Skipping chunk {chunk_index + 1} as no folders exist.")
        continue

    # Generate zip name based on the first and last T in the chunk
    first_T = chunk_TStrs[0]
    last_T = chunk_TStrs[-1]
    zip_name = f"T{first_T}_to_T{last_T}.zip"
    zip_path = os.path.join(csv_T_root, zip_name)

    # Create the zip command
    cmd = ['zip', '-r', zip_name] + existing_folders
    print(f"Creating {zip_name} with folders: {', '.join(existing_folders)}")

    # Execute the command
    try:
        subprocess.run(cmd, cwd=csv_T_root, check=True)
        print(f"Successfully created {zip_name}")
    except subprocess.CalledProcessError as e:
        print(f"Error creating {zip_name}: {e}")