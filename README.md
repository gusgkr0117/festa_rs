# RUST based FESTA

Implementation of the research paper [FESTA: Fast Encryption from Supersingular Torsion Attacks](https://eprint.iacr.org/2023/1747) using Rust.

# Benchmarks
```cmd
$> cargo run --release
=========TEST #0=======
keygen elapsed : 9444 ms
eval elapsed : 3792 ms
inverse elapsed : 6838 ms
=========TEST #1=======
keygen elapsed : 9521 ms
eval elapsed : 3756 ms
inverse elapsed : 6401 ms
=========TEST #2=======
keygen elapsed : 10339 ms
eval elapsed : 4631 ms
inverse elapsed : 6456 ms
=========TEST #3=======
keygen elapsed : 9250 ms
eval elapsed : 3583 ms
inverse elapsed : 6388 ms
=========TEST #4=======
keygen elapsed : 10023 ms
eval elapsed : 3803 ms
inverse elapsed : 6350 ms
=========TEST #5=======
keygen elapsed : 9824 ms
eval elapsed : 3795 ms
inverse elapsed : 6868 ms
=========TEST #6=======
keygen elapsed : 9535 ms
eval elapsed : 3565 ms
inverse elapsed : 6542 ms
=========TEST #7=======
keygen elapsed : 9128 ms
eval elapsed : 3567 ms
inverse elapsed : 6351 ms
=========TEST #8=======
keygen elapsed : 10386 ms
eval elapsed : 4175 ms
inverse elapsed : 6156 ms
=========TEST #9=======
keygen elapsed : 9121 ms
eval elapsed : 3570 ms
inverse elapsed : 6147 ms
```

# References for the implementation

- [FESTA-PKE/FESTA_SageMath](https://github.com/FESTA-PKE/FESTA-SageMath) with its [paper](https://eprint.iacr.org/2023/660)
- [ThetaIsogenies/two-isogenies](https://github.com/ThetaIsogenies/two-isogenies/tree/main) with its [paper](https://eprint.iacr.org/2023/1747)
- [Efficient computation of the image curve](https://eprint.iacr.org/2018/782)
- [GiacomoPope/KummerIsogeny](https://github.com/GiacomoPope/KummerIsogeny/tree/main)
- [pornin/crrl](https://github.com/pornin/crrl)
- [Use of Entangled basis](https://eprint.iacr.org/2017/1143)
