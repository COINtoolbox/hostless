# Author: Lilianne Nakazono
import numpy as np 

def calculate_distance(table_seg:np.ndarray, radius:float=7):
    """Calculates euclidian distance between the transient and the closest masked source (in pixels)

    Args:
        table_seg (np.ndarray): stamp of the segmented image

    Returns:
        int,int,float,bool: index i, index j, euclidian distance (pixels), flag (True if the transient is hostless)
    """
    table_seg[table_seg>0]=1
    transient_idx = 30
    # breaker = False
    mask_i = None
    mask_j = None
    distance = 100 #out-of-range value
    flag = True #if true, the transiet is hostless


    if table_seg[transient_idx][transient_idx]==1:
        return transient_idx, transient_idx, 0, False
    
    # If not, check what is the position of the closest masked source and calculate euclidian distance
    # Note: the first adjacent pixel that is masked will break the for loop (including diagonally adjacent); 
    
    # for step in np.arange(0,30):
    for step in np.arange(0, radius+1):
        array = np.arange(transient_idx-1-step,transient_idx+2+step)
        for indx, i  in enumerate(array):
            if indx==0 or indx==transient_idx+1+step:
                for j in array:
                    if table_seg[i][j]==1:
                        distance = np.sqrt((transient_idx-i)**2+(transient_idx-j)**2)
                        mask_i = i 
                        mask_j = j
                        flag = False
                        # breaker=True
                        return mask_i, mask_j, distance, flag

            else:
                for j in [array[0],array[-1]]:
                    if table_seg[i][j]==1:
                        distance = np.sqrt((transient_idx-i)**2+(transient_idx-j)**2)
                        mask_i = i 
                        mask_j = j
                        
                        # if inside_radius(distance=distance, radius=radius): #check if masked pixel is inside threshold
                        # breaker=True
                        flag = False
                        return mask_i, mask_j, distance, flag

    return mask_i, mask_j, distance, flag