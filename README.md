# Solar_Freshwater_Desalination
 
Hello! This is a MATLAB project which estimates the amount of potable water produced by a solar still over the course of a day. This was a senior design project done in my 4th year of mechanical engineering aimed to design a cheap, reliable and portable freshwater desalination device that can be easily deployed in freshwater-stressed areas.

![freshwater1](https://user-images.githubusercontent.com/80991754/123301203-66143180-d4e9-11eb-8220-4ddf6ce4430f.JPG)

To calculate the amount of freshwater produced by this solar desalination device, ideally a set of three nonlinear differential equations is required. Though by making the assumption that any component of the solar still over a given timestep remains the same temperature (while keeping the timestep very small), this differential equations can be reduced to a set of three nonlinear heat transfer equations which is much easier to solve. The three heat transfer equations can graphically be seen below.

![freshwater2](https://user-images.githubusercontent.com/80991754/123302611-e12a1780-d4ea-11eb-9f2a-bd4f7d1b0d31.JPG)

Solving these equations with given parameters of the region such as solar flux, ambient air temperature, length of day etc... the below plots were generated. Importantly even in worst case conditions the solar still generated 4.15L, enough to sustain the daily drinking water needs of two adults and two children.

![freshwater3](https://user-images.githubusercontent.com/80991754/123302793-0fa7f280-d4eb-11eb-9486-2cdbd0df55c2.JPG)

![freshwater4](https://user-images.githubusercontent.com/80991754/123302802-12a2e300-d4eb-11eb-97ee-54f93bd0aa31.JPG)

Bit more realistic efficiency values when using a lower (realistic) basin absorptivity coefficient...

![freshwater5](https://user-images.githubusercontent.com/80991754/123302814-18002d80-d4eb-11eb-8e8c-7abd13fa6d11.JPG)
