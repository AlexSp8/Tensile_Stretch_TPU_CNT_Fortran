import GEOM
from salome.geom import geomBuilder
import salome
import math

salome.salome_init()
geompy = geomBuilder.New()

# Parameters (mm) for half specimen
Z  = [115.0, 25.0, 33.0]
Z3 = (Z[0]-2*Z[1]-Z[2])/2
Z.append(Z3)

xy_Symmetry=True

X = [3.2]
if xy_Symmetry:
    X[0]=X[0]/2
    
Y = [10.0, 3.0]
R = (Z[3]**2+(Y[0]-Y[1])**2)/(2*(Y[0]-Y[1]))
Y.append(Y[1]+R)

Vertices = []
Lines = []

yi = [0.0, Y[0], Y[0], Y[1], Y[1], Y[0], Y[0], 0.0]
zi = [0.0, 0.0, Z[1], Z[1]+Z[3], Z[1]+Z[3]+Z[2], Z[0]-Z[1], Z[0], Z[0]]
#Straight lines
for i in range(len(yi)):
    vertex=geompy.MakeVertex(0.0, yi[i], zi[i])
    Vertices.append(vertex)
    check = (i != 0) and (i != 3) and (i != 5)
    if check:
        line = geompy.MakeLineTwoPnt(Vertices[i-1], Vertices[i])
        Lines.append(line)

lengthV=len(Vertices)
yi = [0.0, 0.0, 0.0, 0.0]
zi = [Z[0]-Z[1], Z[0]-Z[1]-Z[3], Z[1]+Z[3], Z[1]]
#Straight lines
for i in range(len(yi)):
    vertex=geompy.MakeVertex(0.0, yi[i], zi[i])
    Vertices.append(vertex)
    j=lengthV+i
    line = geompy.MakeLineTwoPnt(Vertices[j-1], Vertices[j])
    Lines.append(line)

line = geompy.MakeLineTwoPnt(Vertices[len(Vertices)-1], Vertices[0])
Lines.append(line)

#Arc lines
yi_arc=[Y[2], Y[2]]
zi_arc=[Z[1]+Z[3], Z[0]-(Z[1]+Z[3])]
#0-indexing
Vertices_arc1=[2,4]
Vertices_arc2=[3,5]
for i in range(len(yi_arc)):
    vertex=geompy.MakeVertex(0.0, yi_arc[i], zi_arc[i])
    Vertices.append(vertex)
    v1=Vertices_arc1[i]
    v2=Vertices_arc2[i]
    line = geompy.MakeArcCenter(vertex, Vertices[v1], Vertices[v2])
    Lines.append(line)

# for i in range(len(Vertices)):
#     geompy.addToStudy(Vertices[i],"Vertices_"+str(i))

# for i in range(len(Lines)):
#     geompy.addToStudy(Lines[i],"Line_"+str(i))

wire = geompy.MakeWire(Lines)
# geompy.addToStudy(wire,"Wire")

face = geompy.MakeFace(wire, isPlanarWanted=True)
# geompy.addToStudy(face,"Face")

#Extrude
direction_vector = geompy.MakeVectorDXDYDZ(1, 0, 0)
solid = geompy.MakePrismVecH(face, direction_vector, X[0])
geompy.addToStudy(solid, "ASTM D638 Type IV")

########################################################################
# make groups in solid
j = 1  #group counter
e = 1e-4

IDs = []
Groups = []

xi = [0.0, X[0]]
#down, down curve, mid, up curve, up
yi = [Y[0], Y[1], Y[1], Y[1], Y[0]]
zi = [+e, Z[1]+Z[3]-e, Z[0]/2, Z[0]-Z[1]-Z[3]+e, Z[0]-e]
#1-10: z
for k in range(len(zi)):
    Group = geompy.CreateGroup(solid, geompy.ShapeType["EDGE"])
    for i in range(len(xi)):

        vertex = geompy.MakeVertex(xi[i], yi[k], zi[k])
        edge = geompy.GetShapesNearPoint(solid, vertex, geompy.ShapeType["EDGE"])
        ID = geompy.GetSubShapeID(solid, edge)
        IDs.append(ID)
        geompy.UnionIDs(Group, [ID])

    geompy.addToStudyInFather(solid, Group, 'Group_'+str(j))
    Groups.append(Group)
    j = j+1

    Group = geompy.CreateGroup(solid, geompy.ShapeType["EDGE"])
    for i in range(len(xi)):

        vertex = geompy.MakeVertex(xi[i], 0.0, zi[k])
        edge = geompy.GetShapesNearPoint(solid, vertex, geompy.ShapeType["EDGE"])
        ID = geompy.GetSubShapeID(solid, edge)
        IDs.append(ID)
        geompy.UnionIDs(Group, [ID])

    geompy.addToStudyInFather(solid, Group, 'Group_'+str(j))
    Groups.append(Group)
    j = j+1

#11: y
Group = geompy.CreateGroup(solid, geompy.ShapeType["EDGE"])
for i in range(len(xi)):
  vertex = geompy.MakeVertex(xi[i], Y[0]/2, 0.0)
  edge = geompy.GetShapesNearPoint(solid, vertex, geompy.ShapeType["EDGE"])
  ID = geompy.GetSubShapeID(solid, edge)
  IDs.append(ID)
  geompy.UnionIDs(Group, [ID])
  vertex = geompy.MakeVertex(xi[i], Y[0]/2, Z[0])
  edge = geompy.GetShapesNearPoint(solid, vertex, geompy.ShapeType["EDGE"])
  ID = geompy.GetSubShapeID(solid, edge)
  IDs.append(ID)
  geompy.UnionIDs(Group, [ID])

geompy.addToStudyInFather(solid, Group, 'Group_'+str(j))
Groups.append(Group)
j = j+1

#12: x
yi = [0.0, Y[0]]
Group = geompy.CreateGroup(solid, geompy.ShapeType["EDGE"])
for i in range(len(yi)):
  vertex = geompy.MakeVertex(X[0]/2, yi[i], 0.0)
  edge = geompy.GetShapesNearPoint(solid, vertex, geompy.ShapeType["EDGE"])
  ID = geompy.GetSubShapeID(solid, edge)
  IDs.append(ID)
  geompy.UnionIDs(Group, [ID])
  vertex = geompy.MakeVertex(X[0]/2, yi[i], Z[0])
  edge = geompy.GetShapesNearPoint(solid, vertex, geompy.ShapeType["EDGE"])
  ID = geompy.GetSubShapeID(solid, edge)
  IDs.append(ID)
  geompy.UnionIDs(Group, [ID])

geompy.addToStudyInFather(solid, Group, 'Group_'+str(j))
Groups.append(Group)
j = j+1

#13: Y=Y0 face
yi = [Y[0], Y[0], Y[1], Y[1], Y[0]]
zi = [0.0, Z[1], Z[1]+Z[3], Z[1]+Z[3]+Z[2], Z[1]+Z[3]+Z[2]+Z[3]]
Group = geompy.CreateGroup(solid, geompy.ShapeType["FACE"])
for i in range(len(zi)):
  vertex = geompy.MakeVertex(X[0]/2, yi[i], zi[i]+e)
  # geompy.addToStudy(vertex,"vertex"+str(i))
  face = geompy.GetShapesNearPoint(solid, vertex, geompy.ShapeType["FACE"])
  ID = geompy.GetSubShapeID(solid, face)
  IDs.append(ID)
  geompy.UnionIDs(Group, [ID])

geompy.addToStudyInFather(solid, Group, 'Group_'+str(j))
Groups.append(Group)
j = j+1


###
### SMESH component
###

import  SMESH, SALOMEDS
from salome.smesh import smeshBuilder

smesh = smeshBuilder.New()

Sub_meshes = []

Mesh_1 = smesh.Mesh(solid,'Mesh_1')

Regular_1D = Mesh_1.Segment()
smesh.SetName(Regular_1D.GetAlgorithm(), 'Regular_1D')

#Default number of segments if it is left unspecified somewhere
Number_of_Segments_1 = Regular_1D.NumberOfSegments(1)
smesh.SetName(Number_of_Segments_1, 'Number_of_Segments_1')

Quadrangle_2D = Mesh_1.Quadrangle(algo=smeshBuilder.QUADRANGLE)
# Quadrangle_2D = Mesh_1.Quadrangle()
smesh.SetName(Quadrangle_2D.GetAlgorithm(), 'Quadrangle_2D')
Quadrangle_Parameters_1 = Quadrangle_2D.QuadrangleParameters(smeshBuilder.QUAD_TRIANGLE_PREF,-1,[],[])
smesh.SetName(Quadrangle_Parameters_1, 'Quadrangle Parameters_1')

# Hexa_3D = Mesh_1.Hexahedron(algo=smeshBuilder.Hexa)
Prism_3D = Mesh_1.Prism()

for i in range(0,len(Groups)):
    Group_ = Mesh_1.GroupOnGeom(Groups[i],'Group_'+str(i+1))
    #smesh.SetName(Group_, 'Group_'+str(i+1))
    # Group_ = Mesh_1.GetGroups()

M_ref = 2.0   #Mesh density

j = 0         #Sub-mesh counter
k = 0         #scale factor

#order: z (down, down curve, mid, up curve, up)
Nelem = [25, 30, 110, 30, 25]
if xy_Symmetry:
    for i in range(len(Nelem)):
        Nelem[i] = 2*Nelem[i]

print('Z down-up: ',Z[1],'mm /',int(M_ref*(Nelem[0])))
print('Z curve: ',Z[3],'mm /',int(M_ref*(Nelem[1])))
print('Z mid: ',Z[2],'mm /',int(M_ref*(Nelem[2])))

Factors = [1.0, 1.0/3.477, 1.0, 3.477, 1.0]
for i in range(len(Factors)):

    Nel = M_ref*(Nelem[i])

    Regular_1D_ = Mesh_1.Segment(geom=Groups[j])
    Number_of_Segments = Regular_1D_.NumberOfSegments(int(Nel))
    smesh.SetName(Number_of_Segments, 'Number of Segments_'+str(j+1))
    Number_of_Segments.SetScaleFactor(Factors[i])
    smesh.SetName(Number_of_Segments, 'Number of Scaled Segments_'+str(j+1))
    # Propagation_of_Node = Regular_1D_.PropagationOfDistribution()
    # smesh.SetName(Propagation_of_Node, 'Propagation of Node Distribution on Opposite Edges')
    j=j+1

    Regular_1D_ = Mesh_1.Segment(geom=Groups[j])
    Number_of_Segments = Regular_1D_.NumberOfSegments(int(Nel))
    smesh.SetName(Number_of_Segments, 'Number of Segments_'+str(j+1))
    Number_of_Segments.SetScaleFactor(1.0/Factors[i])
    smesh.SetName(Number_of_Segments, 'Number of Scaled Segments_'+str(j+1))
    j=j+1

#order: y, x
Nelem = [10, 3]
if xy_Symmetry:
    Nelem = [20, 3]
print('Y: ',Y[1],'-',Y[0],'mm /',int(M_ref*(Nelem[0])))
print('X: ',X[0],'mm /',int(M_ref*(Nelem[1])))

Factors = [1.0, 1.0]
for i in range(len(Factors)):
    Nel = M_ref*(Nelem[i])
    Regular_1D_ = Mesh_1.Segment(geom=Groups[j])
    Number_of_Segments = Regular_1D_.NumberOfSegments(int(Nel))
    smesh.SetName(Number_of_Segments, 'Number of Segments_'+str(j+1))
    Number_of_Segments.SetScaleFactor(Factors[i])
    smesh.SetName(Number_of_Segments, 'Number of Scaled Segments_'+str(j+1))
    # Propagation_of_Node = Regular_1D_.PropagationOfDistribution()
    # smesh.SetName(Propagation_of_Node, 'Propagation of Node Distribution on Opposite Edges')
    j=j+1

########################################################################
#compute the mesh
isDone = Mesh_1.Compute()

Mesh_1.SplitVolumesIntoTetra(Mesh_1, 1)

size_Groups = len(Groups)
face = Mesh_1.GroupOnGeom(Groups[size_Groups-1],'Face')

try:
    Mesh_1.ExportDAT(r'/home/aspyridakis/Downloads/MESH.DTA')
    Mesh_1.ExportDAT(r'/home/aspyridakis/Downloads/FACE.DTA', face, renumber=False)
    # Mesh_1.ExportDAT(r'C:/Users/alex_/Downloads/MESH.DTA')
    pass
except:
    print('ExportPartToDAT() failed. Invalid file name?')

########################################################################
if salome.sg.hasDesktop():
    salome.sg.updateObjBrowser()
