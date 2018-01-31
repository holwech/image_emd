import scipy
import numpy as np
import scipy as sp
import scipy.interpolate
import matlab.engine


# Starts the matlab engine and returns the engine as a object
def start_matlab_engine():
    return matlab.engine.start_matlab()


# Performs the Empirical Mode Decomposition on some bidirectional data.
# + "img" is a numpy array of some size MxN and engine is the matlab engine.
# + "depth" is the sensitivity of the extremas. A higher number results in fewer extremas.
# + "spline_lib": "alglib" OR "scipy". Alglib is recommended as it is significantly
# faster and much more memory efficient.
def IEMD(img, epsilon, max_imfs, engine, depth=0, spline_lib='alglib', rbase=1.0, nlayers=3, lambdaNS=0.0):
    residue = np.copy(img)
    imfs = np.zeros((max_imfs, img.shape[0], img.shape[1]))
    count = 0
    for i in range(0, max_imfs):
        print("Getting IMF ", count + 1)
        imf = sifting(residue, epsilon, engine, depth, spline_lib, rbase, nlayers, lambdaNS)
        residue = residue - imf
        imfs[i,:,:] = imf
        count = count + 1
        print("IMF ", count + 1, " done")
        if monotonic(residue, engine, depth):
            return imfs[:count], residue
    return imfs[:count], residue

# Performs the sifting process until a single IMF is found
def sifting(img, epsilon, engine, depth, spline_lib, rbase, nlayers, lambdaNS):
    h_prev = img
    mean = single_sifting(h_prev, engine, depth, spline_lib, rbase, nlayers, lambdaNS)
    h_curr = h_prev - mean
    count = 0
    while not (sd(h_curr, h_prev) < epsilon):
        mean = single_sifting(h_curr, engine, depth, spline_lib, rbase, nlayers, lambdaNS)
        h_prev = h_curr
        h_curr = h_prev - mean
        count = count + 1
        print("Sifting loop ", count, " sd: ", sd(h_curr, h_prev))
    return h_curr


# Performs a single sifting
def single_sifting(img, engine, depth, spline_lib, rbase, nlayers, lambdaNS):
    maxima, minima, maxima_loc, minima_loc = extrema(img, engine, depth)

    x_max, y_max, z_max = triplex_coords(maxima_loc, maxima)
    x_min, y_min, z_min = triplex_coords(minima_loc, minima)

    upper_img = None
    lower_img = None

    if spline_lib == 'scipy':
        upper_spline = create_spline(x_max, y_max, z_max)
        lower_spline = create_spline(x_min, y_min, z_min)
        _,_, upper_img = reconstruct(img.shape, upper_spline)
        _,_, lower_img = reconstruct(img.shape, lower_spline)
    elif spline_lib == 'alglib':
        upper_spline = create_spline_v2(x_max, y_max, z_max, rbase, nlayers, lambdaNS)
        lower_spline = create_spline_v2(x_min, y_min, z_min, rbase, nlayers, lambdaNS)
        _,_, upper_img = reconstruct_v2(img.shape, upper_spline)
        _,_, lower_img = reconstruct_v2(img.shape, lower_spline)
    else:
        raise Exception("Undefined name ", spline_lib, " for variable spline_lib")

    mean = surface_mean(upper_img, lower_img)

    return mean

# Checks whether the residue is monotonic or not, that is, if the residue has less than 2 extrema.
def monotonic(residue, engine, depth):
    _,_, maxima_loc, minima_loc = extrema(residue, engine, depth)
    num_extrema = np.sum(maxima_loc) + np.sum(minima_loc)
    if num_extrema < 2:
        return True
    return False


# Calculates the stopping criterion value (sd = standard deviation??)
def sd(h_curr, h_prev):
    numerator = np.square(h_prev - h_curr)
    denominator = np.square(h_prev)
    return np.sum(numerator / denominator)

# Calculated the piecewise mean of the upper and lower interpolation
def surface_mean(upper_img, lower_img):
    return (upper_img + lower_img) / 2

# Reconstructs the upper and lower plate based on found spline using SciPy algorithm.
def reconstruct(size, spline):
    x_grid = np.linspace(0, size[0], size[0])
    y_grid = np.linspace(0, size[1], size[1])
    X, Y = np.meshgrid(x_grid, y_grid, indexing='xy')
    Z = np.zeros(size)
    Z = spline(X,Y)
    return X, Y, Z

# Reconstruction for alglib
def reconstruct_v2(size, model):
    x_grid = np.linspace(0, size[0], size[0])
    y_grid = np.linspace(0, size[1], size[1])
    X, Y = np.meshgrid(x_grid, y_grid, indexing='xy')
    Z = xal.rbfgridcalc2v(model, x_grid.tolist(), X.shape[0], x_grid.tolist(), Y.shape[0])
    Z = Z = np.asarray(Z).reshape((75,75))
    return X, Y, Z

def reconstruct_image_old(size, spline):
    new_img = np.zeros(size)
    for iy, row in enumerate(new_img):
        for ix, col in enumerate(row):
            new_img[iy, ix] = spline(ix, iy)
    return new_img



# Converts a 2D array to its (x, y, z) coordinates in 3D space
def triplex_coords(locations, values):
    num_locs = np.sum(locations, dtype=np.int32)
    iterator = 0
    x = np.zeros(num_locs)
    y = np.zeros(num_locs)
    z = np.zeros(num_locs)
    for iy, row in enumerate(locations):
        for ix, col in enumerate(row):
            if col:
                y[iterator] = iy
                x[iterator] = ix
                z[iterator] = values[iy, ix]
                iterator = iterator + 1
    return x, y, z

# Returns the values of the mixima and minima in their respective locations
def extrema(img, engine, depth):
    mat_img = matlab.double(img.tolist())
    maxima_loc = mat2np(engine.imextendedmax(mat_img, depth))
    minima_loc = mat2np(engine.imextendedmin(mat_img, depth))
    maxima = maxima_loc * img
    minima = minima_loc * img
    return maxima, minima, maxima_loc, minima_loc

# A fast conversion from matlab.double array to numpy array
def mat2np(mat):
    return np.array(mat._data).reshape(mat.size, order='F')

# Constructs the spline using the SciPy Rbf algorithm
def create_spline(x,y,z,function='thin-plate'):
    spline = sp.interpolate.Rbf(x,y,z,function=function)
    return spline

# Constructs the spline using the AlgLib RBF Gaussian algorithm
def create_spline_v2(x, y, z, rbase=1.0,nlayers=3,lambdaNS=0.0):
    xy = np.zeros((x.shape[0], 3))
    xy[:,0] = x
    xy[:,1] = y
    xy[:,2] = z
    xy = xy.tolist()
    model = xal.rbfcreate(2, 1)
    xal.rbfsetpoints(model, xy)
    xal.rbfsetalgohierarchical(model, rbase, nlayers, lambdaNS)
    xal.rbfbuildmodel(model)
    return model

def create_spline_old(x,y,z):
    x_grid = np.linspace(0, 75, len(x))
    y_grid = np.linspace(0, 75, len(y))
    X, Y = np.meshgrid(x_grid, y_grid, indexing='xy')
    Z = np.zeros((x.size, z.size))
    spline = sp.interpolate.Rbf(x,y,z,function='thin-plate')
    Z = spline(X,Y)
    return X, Y, Z, spline
