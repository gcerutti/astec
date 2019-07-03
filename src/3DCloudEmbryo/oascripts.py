from time import time
import vtk
import numpy as np

from  array_dict import array_dict
from triangular_mesh import TriangularMesh 



def SetInput(obj, _input):
    if vtk.VTK_MAJOR_VERSION <= 5:
        obj.SetInput(_input)
    else:
        obj.SetInputData(_input)


def image_to_vtk_cell(img,mesh_center=None,coef=1.0,mesh_fineness=1.0):

    start_time = time()
    print "--> Generating vtk mesh from image"

    vtk_mesh = vtk.vtkPolyData()
    vtk_points = vtk.vtkPoints()
    vtk_triangles = vtk.vtkCellArray()
    vtk_cells = vtk.vtkLongArray()
    
    nx, ny, nz = img.shape
    data_string = img.tostring('F')

    reader = vtk.vtkImageImport()
    reader.CopyImportVoidPointer(data_string, len(data_string))
    if img.dtype == np.uint8:
        reader.SetDataScalarTypeToUnsignedChar()
    else:
        reader.SetDataScalarTypeToUnsignedShort()
    reader.SetNumberOfScalarComponents(1)
    reader.SetDataExtent(0, nx - 1, 0, ny - 1, 0, nz - 1)
    reader.SetWholeExtent(0, nx - 1, 0, ny - 1, 0, nz - 1)
    reader.SetDataSpacing(*img.resolution)
    reader.Update()

    considered_cells = np.unique(img)[1:]

    if mesh_center is None:
        mesh_center = np.array(img.resolution)*np.array(img.shape)/2.

    for label in considered_cells:

        cell_start_time = time()

        cell_volume = (img==label).sum()*np.array(img.resolution).prod()
        contour = vtk.vtkDiscreteMarchingCubes()
        SetInput(contour,reader.GetOutput())
        contour.ComputeNormalsOn()
        contour.ComputeGradientsOn()
        contour.SetValue(0,label)
        contour.Update()

        divisions = int(np.ceil(np.power(cell_volume,1/3.)*mesh_fineness))

        decimate = vtk.vtkQuadricClustering()
        SetInput(decimate,contour.GetOutput())
        decimate.SetNumberOfDivisions(divisions,divisions,divisions)
        decimate.SetFeaturePointsAngle(120)
        decimate.Update()

        cell_polydata = decimate.GetOutput()
        
        polydata_points = np.array([cell_polydata.GetPoints().GetPoint(p) for p in xrange(cell_polydata.GetPoints().GetNumberOfPoints())])
        polydata_center = polydata_points.mean(axis=0)
        polydata_points = polydata_center + coef*(polydata_points-polydata_center) - mesh_center

        cell_points = []
        for p in xrange(cell_polydata.GetPoints().GetNumberOfPoints()):
            pid = vtk_points.InsertNextPoint(polydata_points[p])
            cell_points.append(pid)
        cell_points = array_dict(cell_points,np.arange(cell_polydata.GetPoints().GetNumberOfPoints()))

        for t in xrange(cell_polydata.GetNumberOfCells()):
            poly = vtk_triangles.InsertNextCell(3)
            for i in xrange(3):
                pid = cell_polydata.GetCell(t).GetPointIds().GetId(i)
                vtk_triangles.InsertCellPoint(cell_points[pid])
                vtk_cells.InsertValue(poly,label)

        cell_end_time = time()
        print "  --> Cell",label,":",decimate.GetOutput().GetNumberOfCells(),"triangles (",cell_volume," microm3 ) [",cell_end_time-cell_start_time,"s]"

    vtk_mesh.SetPoints(vtk_points)
    vtk_mesh.SetPolys(vtk_triangles)
    vtk_mesh.GetCellData().SetScalars(vtk_cells)

    print "  <-- Cell Mesh      : ",vtk_mesh.GetPoints().GetNumberOfPoints()," Points,",vtk_mesh.GetNumberOfCells()," Triangles, ",len(considered_cells)," Cells"

    end_time = time()
    print "<-- Generating vtk mesh from image      [",end_time-start_time,"s]"

    return vtk_mesh


def vtk_polydata_to_cell_triangular_meshes(polydata):
    mesh = {} 

    polydata_cell_data = polydata.GetCellData().GetArray(0)
    triangle_cell_start_time = time()
    print "  --> Listing triangles"
    print "      - ",polydata.GetNumberOfCells()," triangles"
    polydata_triangles = np.array([[polydata.GetCell(t).GetPointIds().GetId(i) for i in xrange(3)] for t in xrange(polydata.GetNumberOfCells())])   
    triangle_cell_end_time = time()
    print "  <-- Listing triangles            [",triangle_cell_end_time - triangle_cell_start_time,"s]"

    triangle_cell_start_time = time()
    print "  --> Listing triangle cells"
    triangle_cell = np.array([polydata_cell_data.GetTuple(t)[0] for t in xrange(polydata.GetNumberOfCells())],np.uint16)
    triangle_cell_end_time = time()
    print "  <-- Listing triangle cells     [",triangle_cell_end_time - triangle_cell_start_time,"s]"

    start_time = time()
    print "  --> Creating cell meshes "
    for c in np.unique(triangle_cell):
        #print ' --> Cell '+str(c)
        mesh[c] = TriangularMesh()
        cell_triangles = np.arange(polydata.GetNumberOfCells())[np.where(triangle_cell==c)]
        cell_triangle_points = np.array([[polydata.GetCell(t).GetPointIds().GetId(i) for i in xrange(3)] for t in cell_triangles])
        cell_vertices = np.sort(np.unique(cell_triangle_points))


        mesh[c].points = array_dict(np.array([polydata.GetPoints().GetPoint(v) for v in cell_vertices]),cell_vertices).to_dict()
        mesh[c].triangles = array_dict(cell_triangle_points,cell_triangles).to_dict()
        mesh[c].triangle_data = array_dict(np.ones_like(cell_triangles)*c,cell_triangles).to_dict()
    end_time = time()
    print "  <-- Creating cell meshes   [",end_time-start_time,"s]"

    return mesh

def composed_triangular_mesh(triangular_mesh_dict):
    start_time = time()
    print "--> Composing triangular mesh..."

    mesh = TriangularMesh()

    triangle_cell_matching = {}

    mesh_points = np.concatenate([triangular_mesh_dict[c].points.keys() for c in triangular_mesh_dict.keys()])
    mesh_point_positions = np.concatenate([triangular_mesh_dict[c].points.values() for c in triangular_mesh_dict.keys()])
    mesh.points = dict(zip(mesh_points,mesh_point_positions))

    mesh_triangles = np.concatenate([triangular_mesh_dict[c].triangles.values() for c in triangular_mesh_dict.keys()])
    mesh.triangles = dict(zip(np.arange(len(mesh_triangles)),mesh_triangles))

    mesh_cells = np.concatenate([c*np.ones_like(triangular_mesh_dict[c].triangles.keys()) for c in triangular_mesh_dict.keys()])
    triangle_cell_matching = dict(zip(np.arange(len(mesh_triangles)),mesh_cells))

    end_time = time()
    print "<-- Composing triangular mesh     [",end_time-start_time,"]"
    return mesh, triangle_cell_matching



def save_obj_to_cells(mesh, cell_matching, dir_filename):
    
    start_time =time()
    print "--> Saving .obj"

    for c in np.unique(cell_matching.values()):
        objfile = open(dir_filename+"cell_"+str(c)+".obj",'w')
        objfile.write("#\n")

        vertex_index = {}
        for i,v in enumerate(mesh.points.keys()):
            pos = mesh.points[v]
            objfile.write("v "+str(pos[0])+" "+str(pos[1])+" "+str(pos[2])+"\n")
            vertex_index[v] = i
            
        objfile.write("g object"+str(c)+"\n")
        cell_faces = array_dict(cell_matching).keys()[array_dict(cell_matching).values() == c]
        for t in cell_faces:
            tri = mesh.triangles[t]
            objfile.write("f "+str(vertex_index[tri[0]]+1)+"/"+str(c)+" "+str(vertex_index[tri[1]]+1)+"/"+str(c)+" "+str(vertex_index[tri[2]]+1)+"/"+str(c)+"\n")
        objfile.close()
    end_time = time()
    print "<-- Saving .obj [",end_time-start_time,"s]"
    
def loadOBJ(filename):
    verts = []
    face = []

    for line in open(filename, "r"):
        vals = line.replace(',','.').split()
        if len(vals)>0:
            if vals[0] == "v":
                v = map(float, vals[1:4])
                verts.append(v)
            if vals[0] == "f":
                f = map(int, vals[1:4])
                face.append(f)

    return verts, face