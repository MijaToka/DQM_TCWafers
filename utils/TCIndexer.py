from denseIndexer import denseIndexer
import numpy as np

class TCIndexer:
    def __init__(self):
        self.typecodeIndexer_ = {}
        self.maxROC_ = np.empty(0,dtype=int)
        self.maxTrLink_ = np.empty(0,dtype=int) # equivalente a maxErx
        self.maxTCPerLink_ = np.empty(0,dtype=int) # equivalente a maxChPerErx
        self.offsets_ = np.empty(0,dtype=int)
        self.denseIndexers_ = np.empty(0,dtype=denseIndexer)

    def processNewCell(self,typecode,ROC,trLink,trCell):
        if trLink == '-1' or trCell == '-1': return #Skip non connected TCs
        if not typecode in self.typecodeIndexer_:
            self.typecodeIndexer_[typecode]=len(self.typecodeIndexer_)
            self.maxROC_.resize(len(self.typecodeIndexer_))
            self.maxTrLink_.resize(len(self.typecodeIndexer_))
            self.maxTCPerLink_.resize(len(self.typecodeIndexer_))
        
        typecodeidx = self.typecodeIndexer_[typecode]
        
        if typecode[1]=="H": trLink //= 2
            
        self.maxROC_[typecodeidx] = max(self.maxROC_[typecodeidx],ROC+1)
        self.maxTrLink_[typecodeidx] = max(self.maxTrLink_[typecodeidx],trLink+1)
        self.maxTCPerLink_[typecodeidx] = max(self.maxTCPerLink_[typecodeidx],trCell+1)

    def update(self):
        n = len(self.typecodeIndexer_) 
        self.offsets_ = np.zeros(n,dtype=int)
        self.denseIndexers_ = np.array([denseIndexer(3) for _ in range(n)],dtype=denseIndexer)

        for typecodeidx in range(n):
            maxROCs = self.maxROC_[typecodeidx]
            maxLinks = self.maxTrLink_[typecodeidx]
            maxTCPerLink = self.maxTCPerLink_[typecodeidx]
            self.denseIndexers_[typecodeidx].updateRanges([maxROCs,maxLinks,maxTCPerLink])
            if typecodeidx < n-1:
                self.offsets_[typecodeidx + 1] = self.denseIndexers_[typecodeidx].getMaxIndex()
        self.offsets_ = np.cumsum(self.offsets_)

    def getIndexfromTypecode(self,typecode):
        if typecode not in self.typecodeIndexer_:
            raise ValueError("Unable to find typecode in Cell Indexer")
        return self.typecodeIndexer_[typecode]
    
    def getTypecodefronIndex(self,typecodeidx):
        for typecode in self.typecodeIndexer_:
            if self.typecodeIndexer_[typecode]==typecodeidx:
                return typecode
        raise ValueError("Unable to find index that matches a typecode")
    
    def denseIndex(self,typecodeidx,ROC,TrLink,TC):
        #moduleType = self.getTypecodefronIndex(typecodeidx)[1]
        #ROCMult = 8 if moduleType == "H" else 16
        #globalTrLink = ROC*ROCMult + TrLink
        return self.denseIndexers_[typecodeidx].denseIndex([ROC,TrLink,TC])
    
    def maxDenseIdx(self):
        i=len(self.typecodeIndexer_)
        if i == 0:
            return 0
        return self.offsets_[-1] + self.maxTrLink_[-1]*self.maxTCPerLink_[-1]
    
    def getIndexerFor(self,typecodeidx):
        return self.denseIndexers_[typecodeidx]
    
    def getNTLinksfor(self,typecodeidx):
        return self.maxTrLink_[typecodeidx]
    
    def getNTLinksfromTypecode(self,typecode):
        typecodeidx = self.getIndexfromTypecode(typecode)
        return self.getNTLinksfor(typecodeidx)
    
    def getNWordsExpectedfor(self,typecodeidx):
        return self.getIndexerFor(typecodeidx).getMaxIndex()
    
    def getNWordsExpectedfromTypecode(self,typecode):
        typecodeidx = self.getIndexfromTypecode(typecode)
        return self.getNWordsExpectedfor(typecodeidx)    