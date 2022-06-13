<hr>  </hr>

<a href="https://github.com/nickabattista/IB2d"><img src="https://static.wixstatic.com/media/50968c_6e90280106f24ba3ada127d6e1620ea5~mv2.png/v1/fill/w_443,h_319,al_c,q_80,usm_0.66_1.00_0.01/50968c_6e90280106f24ba3ada127d6e1620ea5~mv2.webp" align="right" height="450" width="450" ></a>
<H1> IB2d </H1>

Author: Nicholas A. Battista, Ph.D. <br>
Email: <a href="mailto:battistn[at]tcnj.edu"> battistn[at]tcnj.edu </a> <br>
Website: <a href="http://battistn.pages.tcnj.edu"> http://battistn.pages.tcnj.edu </a> <br>
Department: Mathematics & Statistics (<a href="https://mathstat.tcnj.edu/">TCNJ MATH</a>) <br>
Institution: The College of New Jersey (<a href="https://tcnj.edu/">TCNJ</a>) <br> 

<H4>An easy to use immersed boundary method in 2D, with full implementations in MATLAB and Python that contains over 60 built-in examples, including multiple options for fiber-structure models and advection-diffusion, Boussinesq approximations, and/or artificial forcing. </H4>

<h3 style="color:red;"> If you use this software for research, educational, or recreational purposes, please let Nick Battista (<a href="mailto:battistn[at]tcnj.edu">battistn[at]tcnj.edu</a>) know! </h3>


<hr>  </hr>


<H3>If you use the code for research, please cite the following papers:</H3>

N.A. Battista, A.J. Baird, L.A. Miller, A mathematical model and MATLAB code for muscle-fluid-structure simulations, Integ. Comp. Biol. 55(5):901-911 (2015), <a href="http://www.ncbi.nlm.nih.gov/pubmed/26337187"> LINK </a>

N.A. Battista, W.C. Strickland, L.A. Miller,  IB2d:a Python and MATLAB implementation of the immersed
boundary method,, Bioinspiration and Biomemetics 12(3): 036003 (2017), <a href="http://iopscience.iop.org/article/10.1088/1748-3190/aa5e08/meta"> LINK </a>

N.A. Battista, W.C. Strickland, A. Barrett, L.A. Miller, IB2d Reloaded: a more powerful Python and MATLAB implementation of the immersed boundary method, in press Math. Method. Appl. Sci. 41:8455-8480 (2018) <a href="http://onlinelibrary.wiley.com/doi/10.1002/mma.4708/epdf?author_access_token=HKAwHFmV1yKY6_lY4_I0dU4keas67K9QMdWULTWMo8P3KIzKeMHgO9D_yBVf1ZxhuLjZr3RgM74HKTOZj3MqwU9I9Skl8KVs-2ruPFMgjIXF0QlZful2HU6NM7TQ0wkl"> LINK </a>

<hr>  </hr>

<H3>IB2d Video Tutorials:</H3>

Tutorial 1: <a href="https://youtu.be/PJyQA0vwbgU"> https://youtu.be/PJyQA0vwbgU </a>    
An introduction to the immersed boundary method, fiber models, open source IB software, IB2d​, and some FSI examples!

Tutorial 2:  <a href="https://youtu.be/jSwCKq0v84s"> https://youtu.be/jSwCKq0v84s </a>    
A tour of what comes with the IB2d software, how to download it, what Example subfolders contain and what input files are necessary to run a simulation

Tutorial 3:  <a href="https://youtu.be/I3TLpyEBXfE"> https://youtu.be/I3TLpyEBXfE </a>  
An overview of how to construct immersed boundary geometries and create the input files (.vertex, .spring, etc.) for an IB2d simulation to run using the oscillating rubberband example from Tutorial 2 as a guide.

Tutorial 4: <a href="https://youtu.be/4D4ruXbeCiQ"> https://youtu.be/4D4ruXbeCiQ </a>  
The basics of visualizing data using open source visualization software called <a href="https://wci.llnl.gov/simulation/computer-codes/visit/"> VisIt </a> (by Lawrence Livermore National Labs), visualizing the Lagrangian Points and Eulerian Data (colormaps for scalar data and vector fields for fluid velocity vectors)

<hr> </hr>

<H3> IB2d News</H3>

--> Our software release paper on the recenty added concentration dynamics infrastructure was published in Bioinspiration & Biomimetics. More information can be found at the following link: 

<p style="margin-left:50px; margin-right:50px;">M. Santiago, N.A. Battista, L.A. Miller, S. Khatri, <a href="https://doi.org/10.1088/1748-3190/ac4afa"> Passive concentration dynamics incorporated into the library IB2d, a two-dimensional implementation of the immersed boundary method </a>, Bioinspiration & Biomimetics 17: 036003 (2022) </p>

-->  <a href="https://www.mateasantiago.com"> Matea Santiago </a>  has updated the advection-diffusion solver to a third-order WENO scheme.

--> We have released a semi-automatic meshing tool, <a href="https://github.com/dmsenter89/MeshmerizeMe"> MeshmerizeMe </a>, to help discretize Lagrangian geometries. More information can be found in our software release paper: 

<p style="margin-left:50px; margin-right:50px;">D.M. Senter, D.R. Douglas, W.C. Strickland, S. Thomas, A. Talkington, L.A. Miller, N.A. Battista, <a href="https://doi.org/10.1088/1748-3190/ababb0"> A Semi-Automated Finite Difference Mesh Creation Method for Use with Immersed Boundary Software IB2d and IBAMR </a>, Bioinspiration and Biomimetics 16(1): 016008 (2021)</p>

--> The MATLAB plotting routine is incompatible with MATLAB R2020a and R2020b. If you run into this issue, please consider visualizing the .vtk data with <a href="https://wci.llnl.gov/simulation/computer-codes/visit/"> VisIt </a> or <a href="https://www.paraview.org/"> ParaView </a> instead.   

<hr> </hr>


<H4> VisIt Workaround for macOS 11 (Big Sur) </H4>

-->  Many VisIt users have gotten the following error when trying to open visualization files:

<p align="center"> "The MetaData server running on localhost could not get the file list for the current directory" </p>

--> To circumvent this error, try opening VisIt from the command line. Try the following:

<ol type="1">
  <li> In your terminal, change directory until you're in the directory you want to open (eg - the viz_IB2d folder) </li>
  <li> In your terminal, type: /Applications/VisIt.app/Contents/Resources/bin/visit -gui  (or the correct path based to where you installed VisIt) </li>
</ol>  

<hr> </hr>

