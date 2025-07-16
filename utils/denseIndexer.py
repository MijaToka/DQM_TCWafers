import numpy as np

class denseIndexer:
    def __init__(self, n):
        self.n_ = n
        self.maxIdx_ = 0
        self.maxCoords_ = np.zeros(n)
        #the coords are ordered such that the first is the most impactfull one eg
        # maxNROC, maxNHalf, maxNChannelPerHalf
    
    def updateRanges(self,maxCoords):
        assert self.n_ == len(maxCoords)
        self.maxCoords_ = maxCoords
        acc = 1
        [acc := acc * v for v in self.maxCoords_] 
        self.maxIdx_ = acc

    def denseIndex(self,coords):
        assert len(coords) == len(self.maxCoords_)
        output = coords[0]
        for i in range(1,len(coords)):
            output = output*self.maxCoords_[i] + coords[i]
        return output
    def unpackDenseIndex(self,denseId):
        coords = np.zeros(self.n_)
        for n in range(self.n_-1,-1,-1):
            coords[n] = denseId%self.maxCoords_[n]
            denseId//=self.maxCoords_[n]
        return coords
    def getMaxIndex(self): return self.maxIdx_