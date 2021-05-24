import Util.AMSOS as am


SylinderFileList = am.getFileListSorted("./result*-*/Sylinder_*.pvtp")
ProteinFileList = am.getFileListSorted("./result*-*/Protein_*.pvtp")
ConBlockFileList = am.getFileListSorted("./result*-*/ConBlock_*.pvtp")

print(SylinderFileList)
print(ProteinFileList)
print(ConBlockFileList)

assert len(SylinderFileList) == len(ConBlockFileList)
assert len(SylinderFileList) == len(ProteinFileList)

# example, print 5 files
for f in (SylinderFileList[:5]):
    frame = am.FrameVTK(f)
    frame.printData()

for f in (ProteinFileList[:5]):
    frame = am.FrameVTK(f)
    frame.printData()

for f in (ConBlockFileList[:5]):
    frame = am.FrameVTK(f)
    frame.printData()
