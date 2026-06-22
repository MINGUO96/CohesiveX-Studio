# Cohesive-zone model notes

For a zero-thickness cohesive interface, the displacement jump is

```text
Delta u = u+ - u-
```

The local separations are

```text
delta_n = Delta u dot n
delta_s = Delta u dot s
delta_t = Delta u dot t
delta_sh = sqrt(delta_s^2 + delta_t^2)
```

The initial penalty response is

```text
t_n = K_I delta_n
t_s = K_II delta_s
t_t = K_II delta_t
```

For a bilinear triangular law, the critical separations are

```text
delta_n_f = 2 G_Ic / S_I
delta_s_f = 2 G_IIc / S_II
```

The Benzeggagh-Kenane mixed-mode critical energy is

```text
G_c = G_Ic + (G_IIc - G_Ic) (G_s / G_T)^eta
```

CohesiveX Studio assigns elements to cohesive families using domain relationships:

```text
grain-boundary interface: domain(e1) != domain(e2)
intragranular interface:  domain(e1) == domain(e2)
```

Intragranular sampling uses a deterministic interface-level random number when a seed is provided, or a documented default seed when a fractional selection is requested without an explicit seed.
