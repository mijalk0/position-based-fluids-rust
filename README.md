# position-based-fluids-rust
Position Based Fluids implemented in Rust, method from Müeller and Macklin 2013

# Preview

![](https://github.com/mijalk0/position-based-fluids-rust/blob/master/writeup/images/simulation_end.png)

# Install Rust

Ensure `rustc` and `cargo` are installed. If running Unix/Linux/MacOS, issue

```
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

to install rust if needed. Check the [official rust installation guide](https://www.rust-lang.org/tools/install) to be sure.

# Building

To build the binary, clone the repository and `cd` into it. Then issue 

```
cargo build --release
```

# Running

To run the binary, make sure your current working directory is the repository. The binary needs this to parse the `.obj` files of the spheres. Then simply issue

```
cargo run --release
```

To pass in arguments, add them after the `--release` flag, e.g.

```
cargo run --release n 10000
```

to simulate with 10,000 particles.

# Benchmarking

`cd` into the project directory and issue

```
cargo bench
```

to perform benchmarking with [`criterion`](https://docs.rs/criterion/latest/criterion/).

# Uninstall Rust

To uninstall the rust toolchain, issue 

```
rustup self uninstall
```
