#!python3

import numpy as np
import ntpath  
import argparse
import re
import timeit
# Reading of carto data from 
# https://github.com/stephan1312/SlicerEAMapReader
def TextToFloat( text):
    # First remove any leading and trailing spaces and newlines
    leadingNonDigits = "^[^0-9.-]*"
    trailingNonDigits = "[^0-9.-]*$" 
    text = re.sub(leadingNonDigits, "",text)
    text = re.sub(trailingNonDigits, "",text)
    lines = re.split("\n",text)
    numbers = []
    for i in range(len(lines)):
      lines[i] = re.sub(leadingNonDigits, "",lines[i])
      lines[i] = re.sub(trailingNonDigits, "",lines[i])
      lineSplit = re.split("[ ]+", lines[i])
      thisNumbers = []
      for j in range(len(lineSplit)):
        thisNumbers.append(float(lineSplit[j])) 
      numbers.append(thisNumbers)
    return numbers


def readCartoMesh(filename):
    meshName = ntpath.basename(filename)
    #self.addLog("Reading "+meshName+":")
    section = "none"
    verticesText = ""
    trianglesText = ""
    scalarsText = ""
    attributesText = ""
    scalarLabels = ""
    
    with open(filename, "r", encoding="latin-1") as filehandle:
      for line in filehandle:
        #   if self.abortRequested:
        #     return False
        # Remove trailing newline and trailing and leading spaces
        line = re.sub("[\n]$", "", line)
        line = re.sub("[ ]*$", "", line)
        line = re.sub("^[ ]*", "", line)
        
        if len(line) == 0: # empty line
          continue
        if line[0] == ";": # comment line
          continue  
        if line.find("[GeneralAttributes]") > -1:
          section = "general"
          continue  
        if line.find("[VerticesSection]") > -1:
          section = "vertices"  
          continue  
        if line.find("[TrianglesSection]") > -1:
          section = "triangles" 
          continue  
        if line.find("[VerticesColorsSection]") > -1:
          section = "scalars"
          continue  
        if line.find("[VerticesAttributesSection]") > -1:          
          section = "attributes"          
          continue 
        
        if section == "general":
          # Look for scalar labels
          if line.find("ColorsNames") > -1:
            line = re.sub("^ColorsNames[ ]*=[ ]*", "", line)
            scalarLabels = re.split("[ ]*",line)
        if section == "vertices":
          # remove line number ("0 =")
          line = re.sub("[0-9]*[ ]*=[ ]*", "", line)
          # add "clean" line to string
          verticesText = verticesText+line+'\n'
        if section == "triangles":
          line = re.sub("[0-9]*[ ]*=[ ]*", "", line)
          trianglesText = trianglesText+line+'\n'
        if section == "scalars":
          line = re.sub("[0-9]*[ ]*=[ ]*", "", line)
          scalarsText = scalarsText+line+'\n'
        if section == "attributes":
          line = re.sub("[0-9]*[ ]*=[ ]*", "", line)
          attributesText = attributesText+line+'\n'
    
    verticesLong = TextToFloat(verticesText)
    vertices = []
    vertexnormals = []
    for i in range(len(verticesLong)):
      vertices.append([verticesLong[i][0], verticesLong[i][1], verticesLong[i][2]])
      vertexnormals.append([verticesLong[i][3], verticesLong[i][4], verticesLong[i][5]])
      #if self.abortRequested:
      #  return False
    
    trianglesLong = TextToFloat(trianglesText)
    triangles = []
    for i in range(len(trianglesLong)):
      triangles.append([trianglesLong[i][0], trianglesLong[i][1], trianglesLong[i][2]])
      #if self.abortRequested:
      #  return False
      
    if len(scalarsText) > 0:
      scalars = TextToFloat(scalarsText)
    else:
      scalars = []
    
    if len(attributesText) > 0:
      attributes = TextToFloat(attributesText)   # currently not used
    else:
      attributes = []
    
    #self.addLog("  Read "+str(len(vertices))+" vertices, "+str(len(vertexnormals))+" vertex normals, and "+str(len(triangles))+" triangles.")
    #self.addLog("  Read "+str(len(scalarLabels))+" sets of scalars: "+str(scalarLabels)+".")
    
    meshName = "CARTOmesh_"+re.sub(".mesh$", "", meshName)
    #self.addLog("Creating model "+meshName+":")
    
    #modelNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLModelNode')
    
    return {'v':vertices, 'n':vertexnormals, 't':triangles}
    if not self.CreateMesh(modelNode, vertices, vertexnormals, triangles, scalarLabels, scalars):
      slicer.mrmlScene.RemoveNode(modelNode)
      return False
    
    #self.transformCarto(modelNode)
    
    #modelNode.SetName(meshName)  
    #modelNode.CreateDefaultDisplayNodes()
    
    return True

def readCartoPoints( filename):
    pointsName = ntpath.basename(filename)
    #self.addLog("Reading "+pointsName+":")
    
    #fiducialsNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLMarkupsFiducialNode')
    #fiducialsNode.GetMarkupsDisplayNode().SetVisibility(0)
    points=[]
    data=[]
    with open(filename, "r", encoding="latin-1") as filehandle:
      for line in filehandle:
        #if self.abortRequested:
        #  return False
        # Remove trailing newline and trailing and leading spaces
        line = re.sub("[\n]$", "", line)
        line = re.sub("[ ]*$", "", line)
        line = re.sub("^[ ]*", "", line)
        if len(line) == 0: # empty line
          continue
        lineElements = re.split("[ \t]+", line)
        if lineElements[0] == "VERSION_5_0" or lineElements[0] == "VERSION_4_0":
          pointsName = lineElements[1]
        
        if lineElements[0] == "P":
            pointNr = int(lineElements[2])
            pointX = float(lineElements[4]) 
            pointY = float(lineElements[5]) 
            pointZ = float(lineElements[6])  
            unipolar = float(lineElements[10])
            bipolar = float(lineElements[11])
            lat = float(lineElements[12])
            points.append([pointX,pointY,pointZ])
            data.append([unipolar,bipolar,lat])
            #n = fiducialsNode.AddFiducial(pointX, pointY, pointZ)
            #fiducialsNode.SetNthControlPointLabel(n, "Point # "+str(pointNr)+" in "+pointsName)
            #fiducialsNode.SetNthControlPointDescription(n, "Bipolar "+str(bipolar)+" / Unipolar "+str(unipolar)+" / LAT "+str(lat))
            #fiducialsNode.SetNthControlPointLocked(n, 1)
             
    #self.transformCarto(fiducialsNode)        
    
    #pointsName = "CARTOpoints_"+re.sub("_car.txt$", "", pointsName)
    #fiducialsNode.SetName(pointsName)
    #fiducialsNode.GetMarkupsDisplayNode().SetTextScale(0)
    #fiducialsNode.GetMarkupsDisplayNode().SetVisibility(1)
    
    #self.addLog("Created markup fiducials "+pointsName+".")
    #fiducialsNode.CreateDefaultDisplayNodes()
    return {'p':points,'data':data}
    return True
#https://stackoverflow.com/questions/2924795/fastest-way-to-compute-point-to-triangle-distance-in-3d


def pointsToTriangleV(A,B,C,Points):    
    """
    Projects a numpy array of points to a triangle spand by A, B, C
    A: (1 x 3 ) Point 
    B: (1 x 3 ) Point 
    C: (1 x 3 ) Point 
    Points: (n x 3) Point array
    """
    #Build plane from  normal
    n = np.cross(np.subtract(B,A),np.subtract(C,A))
    n *= 1.0/np.linalg.norm(n)
    AC= np.subtract(C,A)
    AB= np.subtract(B,A)
    dot00 = np.dot(AC,AC)
    dot01 = np.dot(AC,AB)
    dot11 = np.dot(AB,AB)
    d = dot00 * dot11 - dot01 *dot01

    dPlaneV= np.dot(Points,n)-np.dot(A,n)

    t1=np.repeat(n ,len(dPlaneV),axis=0).reshape([3,len(dPlaneV)])  
    pPlaneV = np.add(Points, np.multiply(t1, -dPlaneV ).transpose() )
    ApV= np.subtract(pPlaneV,A)    
    dot02V = np.tensordot(ApV,AC,axes=(1))
    dot12V = np.tensordot(ApV,AB,axes=(1)) 
    uV = (dot11 * dot02V - dot01 * dot12V) /d
    vV = (dot00 * dot12V - dot01 * dot02V) /d
    #Now force the point to be on the on triangle
    uV[uV<0]=0 
    vV[vV<0]=0
    t=uV+vV
    uV[t>1]/=t[t>1]
    vV[t>1]/=t[t>1] 
    # Generate the new point
    newPV=  np.add(A,np.add(np.tensordot(uV ,AC, axes=(0)),   np.tensordot(vV ,AB, axes=(0))))

    tV=np.subtract(Points,newPV)
    distanceV=np.sum(tV*tV,axis=1)
      
    return [newPV, distanceV]

def PointsToSurfaceDistanceV(mesh,points):
    """
    Projects all points in the numpy array points to th mesh
    Mesh: Dictionary containing verties ['v'] and triangles ['t']
    Points: (n x 3) Point array
    """
    numP= np.shape(points)[0]
    mindist=(np.arange(numP,dtype=float)+1)*np.inf
    newPoint=(np.arange(numP*3,dtype=float)+1).reshape([numP,3])
    triangles = np.array(mesh['t'], dtype='int')
    vert = mesh['v']
    for tri in triangles:
        # Loop over all verticies        
        A=vert[tri[0]]
        B=vert[tri[1]]
        C=vert[tri[2]]
        ret= pointsToTriangleV(A,B,C,points)
        t=(mindist>ret[1])
        mindist[t]=ret[1][t]
        newPoint[t,:]=ret[0][t] 
    mindist=np.sqrt(mindist)
    return {'newPoints':newPoint, 'distances':mindist}


def flatten(iterable):
    """
    Flattens the iterable object
    """
    try:
        iterator = iter(iterable)
    except TypeError:
        yield iterable
    else:
        for element in iterator:
            yield from flatten(element)

def Interpolate(sourcePoints,sourceData, targetPoints, ssq,maxDist):
  """
  Interpolate the source data from the source points to target points,
  The weigthed interpolation uses the exponential function exp (-distanc/2*ssq).
  The support of the interpolated is maxDist
  """
  weigths=np.zeros([len(sourcePoints),len(targetPoints)])
  for i in range(len(sourcePoints)):
    tV=np.subtract(targetPoints,sourcePoints[i,:])
    distanceV=np.sum(tV*tV,axis=1)
    [ind,]=np.where(distanceV<maxDist)
    x2=np.sum(tV[ind,:]*tV[ind,:],axis=1)
    weigths[i,ind]=np.exp(-x2/(2*ssq))
  s=np.sum(weigths,axis=0)
  l=weigths/s
  weigths=np.nan_to_num(l)
   
  r=np.dot(sourceData[:,2],weigths)
  r[s==0]=np.nan
  return r




def SaveMeshToVtk(filename,points,triangles,data):
  """
  Crude writer for VTK mesh. prevents dependency of VTK
  """
  f = open(filename, "w")
  f.write("# vtk DataFile Version 2.0 \n")
  f.write("Points\n")
  f.write("ASCII\n")
  f.write("DATASET POLYDATA\n")
  f.write("POINTS {} float\n".format(points.shape[0]))
  for p in points:
    f.write("{} {} {}\n".format(p[0],p[1],p[2])) 
  if (len(triangles)>0):
    f.write("POLYGONS {} {} \n".format(triangles.shape[0],triangles.shape[0]*4))
    for p in triangles:
      f.write("3 {0:d} {1:d} {2:d}\n".format(p[0],p[1],p[2])) 

  if (len(data)==points.shape[0]):
    f.write("POINT_DATA {} \n".format(points.shape[0]))
    f.write("SCALARS sample_scalars float 1\n")
    f.write("LOOKUP_TABLE default\n")
    for p in np.nan_to_num(data,nan=-10000):
      f.write("{}\n".format( p)) 

  f.close()



# If used as main 
if __name__ == '__main__':
  """
  Entery point if used as script
  """
  parser = argparse.ArgumentParser("convert")
  parser.add_argument("Mesh", help="Carto Mesh file")
  parser.add_argument("Car", help="Car file of the mesh ")
  parser.add_argument("--maxDistance",default=7.0, type=float ,help="Maximal distance of points from surface")

  parser.add_argument("Output", help="Vtk file to write the result ")
  args = parser.parse_args()
  fileMesh =args.Mesh #"/app/11-1-1-ReFAMLA.mesh"
  fileCar =args.Car #"/app/11-1-1-ReFAMLA_car.txt"
  filenameOut =args.Output # "/data/11-1-1-ReFAMLA.vtk"
  maxDistance=args.maxDistance
 

  print("Read mesh ...",end='')
  mesh = readCartoMesh(fileMesh)
  print("Done ")
  print("Read Points... ",end='')
  points = readCartoPoints(fileCar)
  print(" Done")
  print("Remove Points with lat -1000 ...", end="")
  t1= points['data']
  pointsData= np.array(points['data'])
  points= np.array(points['p'])
  badPoints = pointsData[:,2]!=-10000
  pointsData=pointsData[badPoints]
  points=points[badPoints]
  print("Done")
  print("Project points to Surface ...", end="")
  ret = PointsToSurfaceDistanceV(mesh,points)
  print("Done")
  print("Interpolate ...",end="")
  s=  maxDistance /(2.0*2.302585)
  inter=Interpolate(ret['newPoints'],pointsData, mesh['v'],s,maxDistance*maxDistance)
  print("Done")
  print("save {} ...".format(filenameOut), end="")
  SaveMeshToVtk(filenameOut,np.array(mesh['v']),np.array(mesh['t'],dtype=int),inter)
  print("Done")
