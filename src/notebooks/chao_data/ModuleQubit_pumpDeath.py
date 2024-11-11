import h5py
import matplotlib.pyplot as plt
import numpy as np
import yaml
from tqdm import tqdm

fileName = "Q2_pumpDeath"
# fileName = "Q3_pumpDeath"
saveDir = "/home/evm9/corral_crowding/src/notebooks/chao_data/"

if __name__ == "__main__":
    data = h5py.File(saveDir + fileName)
    freqList = data["freqList"][()]
    pwrList = data["pwrList"][()]
    glist = data["glist"][()]

    plt.figure(figsize=(14, 3))
    # plt.pcolormesh(freqList/1e9, pwrList, glist.T, shading="auto", cmap="viridis_r", vmin=0.2, vmax=1.1)
    plt.pcolormesh(
        freqList / 1e9,
        pwrList,
        glist.T,
        shading="auto",
        cmap="viridis_r",
        vmin=-0.1,
        vmax=0.9,
    )
    plt.yticks([5, 10, 15, 20, 25])
    plt.xticks([0.5, 1, 1.5, 2])
    plt.colorbar(ticks=[0, 1])
    plt.show()
