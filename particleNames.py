import PartProp.PartPropSvc
from GaudiPython.Bindings import AppMgr 

gaudi = AppMgr()

# get the service from Gaudi:
ppSvc = gaudi.ppSvc() 

# loop over all particle properties:
for pp in ppSvc :  print(pp.particle())
  
# print them as table:
print(ppSvc.all())
