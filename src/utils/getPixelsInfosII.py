from src.utils.getPixelsInfos import getPixelStructure
import os
import numpy as np


def getPixelsInfosII(newPath, configFile, data, pixelsInfos, frameStart):
    if not os.path.exists(newPath):
        os.makedirs(newPath)
    if not os.path.exists(newPath + '/badPixels'):
        os.makedirs(newPath + '/badPixels')
    pixelsInfos = []
    pixelsNRings = np.zeros((1, 20), dtype=int)
    filePixelInfo = open(newPath + "/infoPixels.txt", "w+")
    filePixelInfo.close()
    #get average of coherent pixels:
    allCandidates = [i for i in range(len(pixelsInfos)) if pixelsInfos[i][5] == 1]
    for time in range(frameStart, configFile.FrameEnd):
        candidates = [allCandidates[i] for i in range(len(allCandidates)) if (
                (pixelsInfos[allCandidates[i]][6] - 2 <= time) and (pixelsInfos[allCandidates[i]][6] + 2 >= time))]
        if len(candidates) > 0:
            print(candidates)
