from bmtk.builder.auxi.node_params import positions_list
import numpy as np

def pos_blad():
    """Returns positions_list for CA3b Cells (Basket cell interneurons)
    """
    # Create the possible x,y,z coordinates
    x_start, x_end, x_stride = 0,20,1
    y_start, y_end, y_stride = 0,10,1
    z_start, z_end, z_stride = 1,1.5,2
    x_grid = np.arange(x_start,x_end,x_stride)
    y_grid = np.arange(y_start,y_end,y_stride)
    z_grid = np.arange(z_start,z_end,z_stride)
    xx, yy, zz = np.meshgrid(x_grid, y_grid, z_grid)
    pos_list = np.vstack([xx.ravel(), yy.ravel(), zz.ravel()]).T
    return positions_list(pos_list)

def pos_mpg():
    """Returns positions_list for CA3b Cells (Basket cell interneurons)
    """
    # Create the possible x,y,z coordinates
    x_start, x_end, x_stride = 0,20,1
    y_start, y_end, y_stride = 0,10,1
    z_start, z_end, z_stride = 2,2.5,2
    x_grid = np.arange(x_start,x_end,x_stride)
    y_grid = np.arange(y_start,y_end,y_stride)
    z_grid = np.arange(z_start,z_end,z_stride)
    xx, yy, zz = np.meshgrid(x_grid, y_grid, z_grid)
    pos_list = np.vstack([xx.ravel(), yy.ravel(), zz.ravel()]).T
    return positions_list(pos_list)

