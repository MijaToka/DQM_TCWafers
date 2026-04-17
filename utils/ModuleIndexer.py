
import numpy as np
from warnings import warn as Warn
from denseIndexer import denseIndexer
from FedReadoutSequence import FedReadoutSequence
from utils import resized


class ModuleIndexer:
    def __init__(self):
        self.modFedIndexer_ = denseIndexer(2)
        self.modFedIndexer_.updateRanges([maxCaptureBlockPerFED_:=10,maxECONTperCB_:=12]) # Need to look into this

        self.fedReadoutSequence_ = np.empty(0,dtype=FedReadoutSequence)
        
        self.globalTypesCounter_ = np.empty(0,dtype=int)
        self.globalTypesNTrLinks_ = np.empty(0,dtype=int)
        self.globalTypesNTCs_ = np.empty(0,dtype=int)
        
        self.offsetsModule_ = np.empty(0,dtype=int)
        self.offsetsTrLinks_ = np.empty(0,dtype=int)
        self.offsetsTC_ = np.empty(0,dtype=int)
        self.offsetsData_ = np.empty(0,dtype=int)


        self.maxModuleIdx_ = 0
        self.maxTrLinkIdx_ = 0
        self.maxTCIdx_ = 0

        self.typecodeMap_ = {}

    def reassignTypecodeLocation(self,fedid,cur_modIdx,new_modIdx):
        val = (fedid,cur_modIdx)
        newval = (fedid,new_modIdx)

        for typecode in self.typecodeMap_:
            if self.typecodeMap_[typecode]!=val:
                continue
            self.typecodeMap_[typecode]=newval
            break
        return

    def processNewModule(self,fedid,econtIdx,typecodeidx,nTrLinks,nTCs,typecode):
        if fedid>= (oldsize:=len(self.fedReadoutSequence_)):
            self.fedReadoutSequence_ = resized(self.fedReadoutSequence_,fedid+1,FedReadoutSequence()) # Resize
        
        frs = self.fedReadoutSequence_[fedid]
        frs.id = fedid

        idx = self.modFedIndexer_.denseIndex([0,econtIdx]) # This can be made better
        if idx >= frs.readoutTypes_.size:
            frs.readoutTypes_ = resized(frs.readoutTypes_,idx+1,-1)
        frs.readoutTypes_[idx] = typecodeidx

        if typecodeidx >= self.globalTypesCounter_.size:
            self.globalTypesCounter_.resize(typecodeidx+1)
            self.globalTypesNTrLinks_.resize(typecodeidx+1)
            self.globalTypesNTCs_.resize(typecodeidx+1)
            self.offsetsData_.resize(typecodeidx+1)

        
        self.globalTypesCounter_[typecodeidx] += 1
        self.globalTypesNTrLinks_[typecodeidx] = nTrLinks
        self.globalTypesNTCs_[typecodeidx] = nTCs

        if typecode != '':
            if typecode in self.typecodeMap_:
                Warn(f"Found typecode {typecode} already in the map (fedid,modid) = {self.typecodeMap_[typecode]} overwriting with {(fedid,idx)}") 
            self.typecodeMap_[typecode] = (fedid,idx)

    def finalize(self):
        # Max Indices
        nfeds = self.fedReadoutSequence_.size
        maxModulesIdx = sum(self.globalTypesCounter_)
        maxNTrLinksIdx = np.dot(self.globalTypesCounter_,self.globalTypesNTrLinks_)
        maxNTCIdx = np.dot(self.globalTypesCounter_,self.globalTypesNTCs_)
        # Global offset
        self.offsetsModule_.resize(maxModulesIdx)
        self.offsetsTrLinks_.resize(maxModulesIdx)
        self.offsetsTC_.resize(maxModulesIdx)

        for i in range(1,len(self.globalTypesCounter_)):
            self.offsetsModule_[i] = self.globalTypesCounter_[i-1]
            self.offsetsTrLinks_[i] = self.globalTypesCounter_[i-1] * self.globalTypesNTrLinks_[i-1] 
            self.offsetsTC_[i] = self.globalTypesCounter_[i-1] * self.globalTypesNTCs_[i-1] 

        self.offsetsModule_ = np.cumsum(self.offsetsModule_)
        self.offsetsTrLinks_ = np.cumsum(self.offsetsTrLinks_)
        self.offsetsTC_ = np.cumsum(self.offsetsTC_)

        typecounters = np.zeros(self.globalTypesCounter_.size)

        for fedit in self.fedReadoutSequence_:
            nconn = 0
            fedit.moduleLUT_ = resized(fedit.moduleLUT_,fedit.readoutTypes_.size,-1)
            
            for i in range(len(fedit.readoutTypes_)):
                if fedit.readoutTypes_[i] == -1:
                    continue

                self.reassignTypecodeLocation(fedit.id,i,nconn)
                fedit.moduleLUT_[i]=nconn
                nconn+=1
            fedit.readoutTypes_ = [val for val in fedit.readoutTypes_ if val != -1]

            nmods = len(fedit.readoutTypes_)

            fedit.modOffsets_.resize(nmods)
            fedit.TrLinkOffsets_.resize(nmods)
            fedit.TCDataOffsets_.resize(nmods)
            fedit.enabledTrLink_.resize(nmods)

            for i in range(nmods):
                type_val = fedit.readoutTypes_[i]

                baseMod_offset = self.offsetsModule_[type_val]
                fedit.modOffsets_[i]=baseMod_offset

                baseTrLinkOffset = self.offsetsTrLinks_[type_val]
                internalTrLink_offset = self.globalTypesNTrLinks_[type_val] * typecounters[type_val]
                fedit.TrLinkOffsets_[i] = baseTrLinkOffset + internalTrLink_offset

                baseData_offset = self.offsetsData_[type_val]
                internalData_offset = self.globalTypesNTCs_[type_val] * typecounters[type_val]
                fedit.TCDataOffsets_[i] = baseData_offset + internalData_offset

                fedit.enabledTrLink_[i] = (1 << self.globalTypesNTrLinks_[type_val]) - 1
                typecounters[type_val]+=1
    
    def getIndexandFedForTypecode(self,typecode:str):
        if typecode not in self.typecodeMap_:
            raise ValueError("Did not find typecode in the typecodeMap")
        return self.typecodeMap_[typecode]
                
    def getCellType(self,typecode:str):
        pass
    
    def getIndexForModule(self,fedid,modid):
        return self.fedReadoutSequence_[fedid].modOffsets_[modid]

    def getIndexForEle(self,fedid,econtid):
        modid = self.denseIndexing(fedid,econtid)
        return self.getIndexForModule(fedid,modid)

    def  getIndexForTypecode(self,typecode:str):
        (fedid,modid) = self.getIndexandFedForTypecode(typecode)
        return self.getIndexForModule(fedid,modid)

    def getIndexForModuleLink(self,fedid,modid,trLinkid):
        return self.fedReadoutSequence_[fedid].TrLinkOffsets_[modid] + trLinkid

    def getIndexForEcont(self,fedid,econtid,trLinkid):
        modid = self.denseIndexing(fedid,econtid)
        return self.fedReadoutSequence_[fedid].TrLinkOffsets[modid] + trLinkid
    
    def getIndexForTC(self,fedid,econtid,trLinkid,TC):
        modid = self.denseIndexing(fedid,econtid)
        return self.getIndexForTPGData(fedid,modid,trLinkid,TC)

    def getIndexForTPGData(self,fedid,modid,trLinkid,TCid):
        return self.fedReadoutSequence_[fedid].TCDataOffsets_[modid] + trLinkid * 4 + TCid # The for is 4 TC per TLink

    def getIndexForModuleData(self,typecode:str):
        (fedid,modid) = self.getIndexandFedForTypecode(typecode)
        return self.getIndexForTPGData(fedid,modid,0,0)

    def getNumFEDs(self):
        return sum(1 for fedrs in self.fedReadoutSequence_ if fedrs.readoutTypes_.size != 0)
    
    def getMaxFEDSize(self):
        return self.fedReadoutSequence_.size
    
    def getMaxModuleSize(self):
        return self.maxModuleIdx_
    
    def getNumModules(self,fedid):
        return self.fedReadoutSequence_[fedid].readoutTypes_.size
