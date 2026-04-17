import numpy as np

class FedReadoutSequence:
    def __init__(self):
        self.id : int = 0
        self.moduleLUT_ = np.empty(0,dtype=int)
        self.readoutTypes_ = np.empty(0,dtype=int)
        self.modOffsets_ = np.empty(0,dtype=int)
        self.TrLinkOffsets_ = np.empty(0,dtype=int)
        self.TCDataOffsets_ = np.empty(0,dtype=int)
        self.enabledTrLink_ = np.empty(0,dtype=int)
    
    def __str__(self):
        return f"FedReadoutSequence(id={self.id}, moduleLUT_={self.moduleLUT_}, readoutTypes={self.readoutTypes}, \
modOffsets_={self.modOffsets_}, TrLinkOffsets_={self.TrLinkOffsets_}, TCDataOffsets_={self.TCDataOffsets_}, enabledErx_={self.enabledErx_})"
    
    def __repr__(self):
        return f"FedReadoutSequence(id={self.id}, moduleLUT_={self.moduleLUT_}, readoutTypes={self.readoutTypes}, \
modOffsets_={self.modOffsets_}, TrLinkOffsets_={self.TrLinkOffsets_}, TCDataOffsets_={self.TCDataOffsets_}, enabledErx_={self.enabledErx_})"