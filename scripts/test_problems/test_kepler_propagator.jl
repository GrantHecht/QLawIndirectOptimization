
using QLawIndirectOptimization
using StaticArrays

function main()
    # Define initial and target orbital elements
    r0 = SA[1131.340, -2282.343, 6672.423]
    v0 = SA[-5.64305, 4.30333, 2.42879]

    # Solve
    QLawIndirectOptimization.kepler(r0, v0, 40*60.0, 398600.4418)
end

main()