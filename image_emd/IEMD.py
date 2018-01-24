# Starts the matlab engine and returns the engine as a object
def start_matlab_engine():
    return matlab.engine.start_matlab()


# Performs the Empirical Mode Decomposition on some bidirectional data.
# img is a numpy array of some size MxN and engine is the matlab engine.
# depth is the sensitivity of the extremas. A higher number results in fewer extremas.
def IEMD(img, engine, depth=0):
    



def sifting(img, engine, depth=0):
    

def single_sifting(img, engine, depth=0):
    maxima, minima, maxima_loc, minima_loc = extrema(img, engine)

    x_max, y_max, z_max = 3D_coords(maxima, maxima_loc)
    x_min, y_min, z_min = 3D_coords(minima, minima_loc)

    _,_,_, upper_spline = create_spline(x_max, y_max, z_max)
    _,_,_, lower_spline = create_spline(x_min, y_min, z_min)

    upper_img = reconstruct_image(img.size, upper_spline)
    lower_img = reconstruct_image(img.size, lower_spline)

    mean = surface_mean(upper_img, lower_img)

    return mean

# Checks whether the stopping criterion for and IMF is satisfied
def stopping_criterion():


# Calculated the piecewise mean of the upper and lower interpolation
def surface_mean(upper_img, lower_img)
    return (upper_img + lower_img) / 2

# Reconstructs the upper and lower plate based on found spline.
def reconstruct_image(size, spline):
    new_img = np.zeros(size)
    for iy, row in enumerate(new_img):
        for ix, col in enumerate(row):
            new_img[iy, ix] = spline(ix, iy)
    return new_img


# Converts a 2D array to its (x, y, z) coordinates in 3D space
def 3D_coords(locations, values):
    x = np.array([])
    y = np.array([])
    z = np.array([])
    for iy, row in enumerate(locations):
        for ix, col in enumerate(row):
            if col:
                y = np.append(y, iy)
                x = np.append(x, ix)
                z = np.append(z, values[iy, ix])
    return x, y, z

# Returns the values of the mixima and minima in their respective locations
def extrema(img, engine, depth=0):
    mat_img = matlab.double(img)
    maxima_loc = mat2np(eng.imextendedmax(img, depth))
    minima_loc = mat2np(eng.imextendedmin(img, depth))
    maxima = maxima_loc * img
    minima = minima * img
    return maxima, minima, maxima_loc, minima_loc

# A fast conversion from matlab.double array to numpy array
def mat2np(mat):
    return np.array(mat._data).reshape(mat.size, order='F')

def create_spline(x,y,z):
    x_grid = np.linspace(0, 75, len(x))
    y_grid = np.linspace(0, 75, len(y))
    X, Y = np.meshgrid(x_grid, y_grid, indexing='xy')
    Z = np.zeros((x.size, z.size))
    spline = sp.interpolate.Rbf(x,y,z,function='thin-plate')
    Z = spline(X,Y)
    return X, Y, Z, spline
