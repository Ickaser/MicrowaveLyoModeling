# Average E field, glass, 9GHz
glass_E = Table([
(angle=15u"°", E1=13021.0u"V/m", E2=11070.0u"V/m", E3=12188.0u"V/m"),
(angle=45u"°", E1=10554.0u"V/m", E2=12913.0u"V/m", E3=15114.0u"V/m"),
(angle=75u"°", E1=13288.0u"V/m", E2=10076.0u"V/m", E3=14265.0u"V/m"),
])
# Total power absorbed, glass, 9GHz
glass_P = Table([
(angle=15u"°", P1=4.8782u"W", P2=3.6805u"W", P3=4.6771u"W"),
(angle=45u"°", P1=3.1914u"W", P2=4.9988u"W", P3=6.3824u"W"),
(angle=75u"°", P1=5.1490u"W", P2=2.8527u"W", P3=6.0795u"W"),
])
# Average E field, ice, 9GHz
ice_E = Table([
(angle=15u"°", E1=16041.0u"V/m", E2=12051.0u"V/m", E3=16298.0u"V/m"),
(angle=45u"°", E1=12955.0u"V/m", E2=14115.0u"V/m", E3=16845.0u"V/m"),
(angle=75u"°", E1=12352.0u"V/m", E2=12620.0u"V/m", E3=12944.0u"V/m"),
])
# Total power absorbed, ice, 9GHz
ice_P = Table([
(angle=15u"°", P1=0.17831u"W", P2=0.10242u"W", P3=0.18349u"W"),
(angle=45u"°", P1=0.12012u"W", P2=0.14421u"W", P3=0.20842u"W"),
(angle=75u"°", P1=0.11227u"W", P2=0.11147u"W", P3=0.11729u"W"),
])

glass_E_avg = sum(glass_E.E1 +glass_E.E2 + glass_E.E3) / (3*length(glass_E))
#1.2499e4u"V/m"
ice_E_avg = sum(ice_E.E1 +ice_E.E2 + ice_E.E3) / (3*length(ice_E))
#1.4024e4u"V/m"
glass_P_avg = sum(glass_P.P1 +glass_P.P2 + glass_P.P3) / (3*length(glass_P))
#4.5644u"W"
ice_P_avg = sum(ice_P.P1 +ice_P.P2 + ice_P.P3) / (3*length(ice_P))
#0.142u"W"

total_P = 40u"W"/3

Bf = ice_E_avg^2 / total_P |> u"Ω/m^2"
#1.475e7u"Ω/m^2"
Bvw = glass_E_avg^2 / total_P |> u"Ω/m^2"
#1.17e7u"Ω/m^2"
@show Bf Bvw 