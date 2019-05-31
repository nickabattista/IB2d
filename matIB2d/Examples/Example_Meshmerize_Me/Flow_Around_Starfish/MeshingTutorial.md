# Meshing Tutorial - From Image
## Finding a Good Image
Find a nice image, like this starfish:

![](https://upload.wikimedia.org/wikipedia/commons/f/fd/Reef0297.jpg)

I found it on WikiMedia Commons.

## Getting the Initial Outline

We next start ContourizeMe on the image:

```bash
~$ ContourizeMe ~/sample_images/Reef0297.jpg 
```
Our goal here is to now play with the levers until the contour seems to match our starfish pretty well. (Played with slot on the left-hand side. Need explanation of how the various buttons work.

## Cleaning Up

The next step is to open the SVG created by ContourizeMe with a vector-graphics program such as [Illustrator](https://www.adobe.com/products/illustrator.html) or [Inkscape](https://inkscape.org/). Once opened, we want to look for issues 

![](starfish_initial.svg)

While the initial image is certainly a very close match to the original image, 


## Verifying the Mesh

![](starfish_mesh_image.png)