# pthreads-diffusion-solver
A parallel implementation of the [Finite Difference Method (FDM)](https://en.wikipedia.org/wiki/Finite_difference_method) for solving the 2D heat problem, made using [Pthreads](https://en.wikipedia.org/wiki/Pthreads) in [C](https://en.wikipedia.org/wiki/C_(programming_language)).

![The diffusion solver in action]()

## Setup
### Downloading the repository
```sh
$ git clone X
$ cd X
```

### Installing dependencies
**gnuplot**

Linux/Ubuntu:
```sh
$ sudo apt update
$ sudo apt install gnuplot
```

MacOSX:
```sh
$ brew update
$ brew install gnuplot
```

**ffmpeg**

Linux/Ubuntu:

```sh
$ sudo apt update
$ sudo apt install ffmpeg
```

MacOSX:

```sh
$ brew update
$ brew install ffmpeg
```

### Setting up folder structure
```sh
$ make setup 
```
Creates the folders `data`, `plots`, and `video` which are necessary to run the program. These can be emptied with `make purge`.

## Run
```sh
$ make diffusion_solver
```
Compiles the program. To delete it, use `make clean`.

```sh
$ ./diffusion_solver
```
Runs the program with default arguments.

```sh
$ make video
```
Creates `plots` for the previous run and a `video` animation of them, which is located at `./video/animation.mp4` by default.

## Options
`./diffusion_solver -y [y_size] -x [x_size] -i [iterations] -s [snapshot_freq]`
Option | Description | Restrictions | Default value
:------------ | :------------ | :------------ | :------------
**-y** | Height in pixels | > 0 | 256
**-x** | Width in pixels | > 0 | 256
**-i** | Number of iterations | > 0 | 100000
**-s** | Number of iterations between each snapshot | > 0 | 1000