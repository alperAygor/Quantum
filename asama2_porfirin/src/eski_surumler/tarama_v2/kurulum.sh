#!/bin/bash
# Bulut makinesinde bir kez calistir: bash kurulum.sh
set -e
sudo apt-get update -y && sudo apt-get install -y python3-venv python3-pip
python3 -m venv ~/qenv
source ~/qenv/bin/activate
pip install --upgrade pip
pip install numpy scipy pyscf pennylane tblite rdkit ase h5py
echo "Kurulum tamam. Calistirmak icin: source ~/qenv/bin/activate && bash calistir.sh 4"
