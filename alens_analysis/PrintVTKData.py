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

<<<<<<< HEAD:alens_analysis/PrintVTKData.py
class ConBlock(object):
    end0 = None
    end1 = None
    pass


class Frame:

    def __init__(self, sylinderFile=None, proteinFile=None,
                 conBlockFile=None, print_flag=True):
        self.print_flag = print_flag
        self.sylinders = []
        self.proteins = []
        self.conBlocks = []
        self.proteinBlocks = []
        # self.parseSylinderFile(sylinderFile)
        # self.parseProteinFile(proteinFile)
        self.parseConBlockFile(conBlockFile)

    def parseFile(self, dataFile, objType, objList):
        print("Parsing data from " + dataFile)
        # create vtk reader
        reader = vtk.vtkXMLPPolyDataReader()
        reader.SetFileName(dataFile)
        reader.Update()
        data = reader.GetOutput()

        # fill data
        # step 1, end coordinates
        nObj = int(data.GetPoints().GetNumberOfPoints() / 2)
        print("parsing data for ", nObj, " object(s)")
        for i in range(nObj):
            s = objType()
            s.end0 = data.GetPoints().GetPoint(2 * i)
            s.end1 = data.GetPoints().GetPoint(2 * i + 1)
            objList.append(s)

        # step 2, member cell data
        numCellData = data.GetCellData().GetNumberOfArrays()
        print("Number of CellDataArrays: ", numCellData)
        for i in range(numCellData):
            cdata = data.GetCellData().GetArray(i)
            dataName = cdata.GetName()
            # print("Parsing Cell Data", dataName)
            for j in range(len(objList)):
                setattr(objList[j], dataName, cdata.GetTuple(j))

        # step 3, member point data
        numPointData = data.GetPointData().GetNumberOfArrays()
        print("Number of PointDataArrays: ", numPointData)
        for i in range(numPointData):
            pdata = data.GetPointData().GetArray(i)
            dataName = pdata.GetName()
            # print("Parsing Point Data", dataName)
            for j in range(len(objList)):
                setattr(objList[j], dataName + "0", pdata.GetTuple(2 * j))
                setattr(objList[j], dataName + "1", pdata.GetTuple(2 * j + 1))

        # output all data for debug
        if self.print_flag:
            for s in objList[:10]:
                # print(s.end0, s.end1)
                attrs = vars(s)
                print('*************************************')
                print('\n'.join("%s: %s" % item for item in attrs.items()))
                print('*************************************')

        print("-------------------------------------")

    def parseSylinderFile(self, sylinderFile):
        self.parseFile(sylinderFile, Sylinder, self.sylinders)

    def parseProteinFile(self, proteinFile):
        self.parseFile(proteinFile, Protein, self.proteins)

    def parseConBlockFile(self, conBlockFile):
        self.parseFile(conBlockFile, ConBlock, self.conBlocks)


def getFrameNumber_lambda(filename): return int(
    re.search('_([^_.]+)(?:\.[^_]*)?$', filename).group(1))


def doSomethingPerFrame(frame):
    pass


def main():
    # get file list
    SylinderFileList = glob.glob(
        './result/result*/Sylinder_*.pvtp')
    ProteinFileList = glob.glob(
        './result/result*/Protein_*.pvtp')
    ConBlockFileList = glob.glob(
        './result/result*/ConBlock_*.pvtp')

    # sort as numerical order
    SylinderFileList.sort(key=getFrameNumber_lambda)
    ProteinFileList.sort(key=getFrameNumber_lambda)
    ConBlockFileList.sort(key=getFrameNumber_lambda)

    print(SylinderFileList)
    print(ProteinFileList)
    print(ConBlockFileList)

    assert len(SylinderFileList) == len(ConBlockFileList)
    assert len(SylinderFileList) == len(ProteinFileList)

    # example
    for i in range(len(SylinderFileList)):
        # get frame
        frame = Frame(SylinderFileList[i], ProteinFileList[i],
                      ConBlockFileList[i]
                      )
        # do something
        doSomethingPerFrame(frame)


if __name__ == '__main__':
    main()
=======
for f in (ConBlockFileList[:5]):
    frame = am.FrameVTK(f)
    frame.printData()
>>>>>>> main:amsos_analysis/PrintVTKData.py
