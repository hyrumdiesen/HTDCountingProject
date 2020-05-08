#!/usr/bin/env python
# coding: utf-8

# In[2]:


import numpy as np
import matplotlib.pyplot as plt
import scipy.misc
import scipy
import scipy.ndimage as npimage
import imageio
import re
import math


# In[3]:


# this was the best solution I could find for finding the coordinates from the rectangle data
'''this is a bit of a hack but we convert each element in the list to a string and seach for integers using the re module'''

def find_centers(rects): 
    list_coor = []
    counter = 0
    for i in rects:
        
        r = i[0]
        c = i[1]
        strr = str(r)
        strc = str(c)



        a = strr.find(',') #finds commas where numbers will be
        b = strc.find(',')
        strr2 = strr[a:]
        strc2 = strc[b:]
        intr1 = int(re.search(r'\d+', strr).group()) #each of these finds integers in a string
        intr2 = int(re.search(r'\d+', strr2).group())

        intc1 = int(re.search(r'\d+', strc).group())
        intc2 = int(re.search(r'\d+', strc2).group())

        center_r = (intr1+intr2)/2
        center_c = (intc1+intc2)/2

        center_coordinates = center_r,center_c #converts recrtangle coordinates to the coordinates of the center of the rectangles
        list_coor.append(center_coordinates)
        counter +=1
    return list_coor
    


# In[4]:


def find_dist(p1,p2):
    '''simple formula to find the distance between two coordinates, we will need this later'''
    distance = math.sqrt( ((p1[0]-p2[0])**2)+((p1[1]-p2[1])**2) )

    return(distance)


# In[5]:


def find_all_dist(list1,list_n,maxd=-1):
    '''given list1 of points, will find the closest point in list2 for that particular point, and gives distance'''
    #NUCLEI COORDINATES MUST BE SECOND LIST, NOT FIRST
    short_dist_info_total = []
    for i in list1:
        dist_list = []
        p1 = i
        for j in list_n: #for loop within a for loop, computationally expensive and could be optimized, but it works and is the easiest solution to find closest nuclei 
            p2 = j
            d = find_dist(p1,p2)
                    
            dist_list.append(d)
        smallest = min(dist_list) #finds smallest distance nuclei in the list
        p = 0 #counts number of elements in list to find index of smallest distance nuclei
        for k in dist_list:
            
            if k == smallest:
                dist_index = list_n[p]
                
            p+=1
          
        short_dist_info = [i,dist_index,smallest]
        if maxd != -1: #if distance is not listed or user specifies -1, it ignores checking max distance
            if smallest <=maxd: # checks if distance is smaller than optional max distance maxd argument, if it's bigger, it's not added to the list
                short_dist_info_total.append(short_dist_info) 
        else:
            short_dist_info_total.append(short_dist_info) 
    return short_dist_info_total# [list1coordinate, list2coordinate, distance] # returns the coordinate of the point that a point on list1 was closest with on list2, as well as the distance this was


# In[6]:


def list_sort(r_list,g_list): # will need this later
    '''if value from r_list1 is found on g_list, removes both matching values from lists and adds them to a list of both r and g'''
    g_r_list = []
    for i in r_list:
        for j in g_list: #for loop within a for loop to search for matching values in lists
            if i ==j:
                g_r_list.append(j)
                r_list.remove(i) #important, if a nucleus is found to be both red and green, it's coordinates are removed from each red and green list
                g_list.remove(j)
    return r_list,g_list,g_r_list     
    


# In[7]:


def rect_obj(image, filt):
    '''calls on functions created above to process an image to find the center coordinates of all blobs defined by enclosing rectangles'''
    blobs_filt = (image>filt).astype(int)
    
    labels = scipy.ndimage.label(blobs_filt)[0]
    rect = scipy.ndimage.find_objects(labels)
    centers = find_centers(rect)
    return centers


# In[8]:


def assign_marks(dist_info): #will need this later
    '''takes distance information from the find_all_dist function and adds the nuclei coordinates to their own list for that color'''
    colored_nuclei = []
    for i in dist_info:
        nuc_coor = i[1]
        colored_nuclei.append(nuc_coor)
    return colored_nuclei


# In[9]:


def cat_nuclei(image_g, g_filt, image_n, n_filt, maxd = -1):
    '''uses functions above to take two images with filter and max distance values and add matching nuclei to a color list'''
    #NUCLEI IMAGE MUST BE SECOND, NOT FIRST
    g_centers = rect_obj(image_g,g_filt)
    n_centers = rect_obj(image_n,n_filt)
    
    g_nuclei = find_all_dist(g_centers,n_centers,maxd)
    g_nuclei = assign_marks(g_nuclei)
    
    g_nuclei = list(set(g_nuclei)) #very important, removes recurring nuclei coordinates from the list
    return g_nuclei


# In[ ]:


def final(image_g, g_filt, image_r, r_filt,image_n, n_filt, maxd = -1):
    #putting everything together
    num_nuclei = 0
    rect = rect_obj(image_n, n_filt)
    for i in rect: #counts total number of nuclei using the rect_obj function and counting number of elements
        num_nuclei +=1
        
    g_nuclei = cat_nuclei(image_g, g_filt, image_n, n_filt, maxd) #referencing cat_nuclei function
    r_nuclei = cat_nuclei(image_r, r_filt, image_n, n_filt, maxd)
    
    sorted_cats = list_sort(r_nuclei,g_nuclei)
    num_g_nuclei = 0
    num_r_nuclei = 0
    num_r_g_nuclei = 0
    for i in sorted_cats[0]: #counters to count the number of nuclei in final lists
        num_g_nuclei +=1
    
    for i in sorted_cats[1]:
        num_r_nuclei +=1
    
    for i in sorted_cats[2]:
        num_r_g_nuclei +=1
    total_stained = num_g_nuclei + num_r_nuclei + num_r_g_nuclei
    print('TOTAL CELLS:', num_nuclei) #distplays information of number of cells within each category
    print('Green Cells:', num_g_nuclei)
    print('Red Cells:', num_r_nuclei)
    print('Green and Red Cells:', num_r_g_nuclei)
    print('None Cells:', num_nuclei - total_stained) #counts unstained nuclei by subtracting total cells from total stained cells
    
    File = open("cell_count.txt",'w')#writes file containing to same directory as python file to store information
    File.write("Number of Total Cells: ")
    File.write(str(num_nuclei)+"\n")
    File.write("Number of Green Cells: ")
    File.write(str(num_g_nuclei)+"\n")
    File.write("Number of Red Cells: ") 
    File.write(str(num_r_nuclei)+"\n")
    File.write("Number of Green and Red Cells: ")
    File.write(str(num_r_g_nuclei) + "\n")
    File.write("Number of None Cells: ")
    File.write(str(num_nuclei - total_stained))
    File.close() 


# In[ ]:


g_name = input('Enter the file name of the green channel tif file:')


# In[ ]:


r_name = input('Enter the file name of the red channel tif file:')


# In[ ]:


n_name = input('Enter the file name of the blue (nuclei) channel tif file:')


# In[5]:


img_test_n = plt.imread(n_name) #blue channel of image taken in my lab, exported directly from the microscope and cropped


# In[7]:


#img_test_n_filt = (img_test_n>90).astype(int) #filtering to a value I'm most comfortable with the program using
 #90 seems like a good value here


# In[8]:


img_test_g = plt.imread(g_name) #green channel


# In[9]:


#img_test_g_filt = (img_test_g>80).astype(int)
# 80 seems like a good value here


# In[10]:


img_test_r = plt.imread(r_name) # red channel, much fewer discrete boundaries than the other examples


# In[11]:


#110 might be a good value, but it's hard to say


# In[12]:


g_filt = input('Please enter filter number for green channel:')
g_filt = float(g_filt)


# In[13]:


r_filt = input('Please enter filter number for red channel:')
r_filt = float(r_filt)


# In[ ]:


n_filt = input('Please enter filter number for nuclei (blue) channel:')
n_filt = float(n_filt)


# In[ ]:


dmax = input('(Optional) enter max distance (in number of pixels) will count as a cell marker:')


# In[ ]:


final(img_test_g, g_filt, img_test_r, r_filt, img_test_n, n_filt)


# In[ ]:




