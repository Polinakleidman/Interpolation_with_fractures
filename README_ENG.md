# Interpolation of oil well parameters considering Earth's crust faults.
## Description:
It’s known that many parameters of Earth’s crust  (such as depth of layers, pressure, porosity) in different points can be estimated by values in certain drilled wells. Dependence is build on shortest ways to drilled ways. The main problem of this search is that such dependence disappear if the way intersects the crust fault. In order to fix it, we have to find shortest way enveloping faults.  

## Problem to solve:
### Input:
* rectangular grid of the plane, values of parameters should be found in all nodes
* faulst description (arrays of broken lines)
* coordinates and parameters of drilled wells (they do not have to locate in grid nodes)
### Output:
* values of perameters in every grid node

## Solution:
 *Notice: Considering that grid is dense enough, it is acceptable to use Manhattan or Chebyshev metric instead of Euclid metric. (in my solution I used Chebyshev metric, (Manhattan metric provides 2 times less edges but works worse, later there will be clarification diagrams)*
* To begin with let's describe graph:
  * xleft, ybottom – coordinates of left bottom point of a grid
  * dx, dy - length and width of one cell
  * nx,ny – count of nodes by the OX and OY (total count of nodes: nx*ny, width of grid=(nx – 1)*dx, length=(ny – 1)*dy
  * Data – array of drilled wells, count_data = their amount
  * Faults – faults(array of broken lines) </br>

So, in our graph **nodes** are gril nodes and drilled wells (их nx * ny + count_data) 
* If the node is a well, there ar maximum 4 edges into 4 grid points surrounding it (in case a fault intersects new edge, it is excluded from graph)
* If the node is a grid node, edges to *wells* in square [x -dx, d+dx] * [y -dy, d+dy] are already added. We add maximum 8 more edges to nearest grid nodes (like in the previous case if these edges intersect faults, they are excluded.

  *Let's notice that there are no edges between wells  on purpose, as parameters of wells in one cell can be on average)*

Now the graph is build, every **edge** contains its length, let's start Djikstra from every well now.

**After every Djikstra start checking existance of straight way not intersecting faults from every grid point into current well.(in case of it's existance, taking it as a miminun distance) then it will be shown, why this is important! let's call this updating *BEAUTIFUL***

#### Let's count asymptotics:
Graph building works for O(n+m) = (nx * ny+ count_data + 8 * (nx * ny+count_data) =O(nx * ny+count_data)
Djikstra  O(mlogn) = O(nlogn)
As alomost always nx * ny < 1000 * 1000, count_data < 1000, => total time < O(10^9 * log(10^9)) ~ several minutes

### Let's fo in interpolation.
*Example parameter - depth of upper border of appropriate layer*
* Let's have a look on isolines 
* **Photo 1**
  ![Photo 1](https://github.com/Polinakleidman/Interpolation_with_fractures/blob/main/1.jpg "")
  * 1ts photo – how would the map look line, if the faults didn't exist (faults are red) It's seen, that all isolines are smooth. Current aim is to save smoothness for cases with faults.
  
  * Second picture. It should be noticed that upper graphic contains more *smmoth* isolines than lower one, but its isolines are diamond-shaped, which against nature lows, at the same time lowest graphic includes a lot of *cloves*, which should be removed. These 2 photos' code differ only in including  ***BEAUTIFUL*** update.
  * **Photo 2**
  ![Photo 2](https://github.com/Polinakleidman/Interpolation_with_fractures/blob/main/no_diag.jpg "") 
  
  * Third picture. It's analogy to the second one, but now Chebyshev metric is used instead of Manhattan. As we see, lowest graphic is almost perfect: *cloves* are gone, isolines are round-shaped, but there are still *vibrations*.
  * **Photo 2b**
  ![Photo 3](https://github.com/Polinakleidman/Interpolation_with_fractures/blob/main/diag.jpg "")
  * How to get rid of vibrations:
    
    In formula:
    ```
      h = sqrt(dxn * dxn + dyn * dyn + asmoosnes * asmoosnes);
      ph = 1 / pow(h , power); 
    ``` 
      numbers for *power* can lead to full smoothing. 
