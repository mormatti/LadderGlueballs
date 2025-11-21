for g in ["ZZ3", "SU3"]
    println("Group: ", g)
    for L in [11]
        println("   Length: ", L)
        for λ in [1//10, 5//10, 9//10]
            println("       Lambda: ", λ)
            for B in ["1"]
                println("           Band: ", B)
                for ℓ in [0,1,2]
                    println("               Support half: ", ℓ)
                    prec = creator[g][λ][L][B][ℓ]["precision"]
                    preceff = creator[g][λ][L][B][ℓ]["precisioneff"]
                    println("                   Precision: ", prec)
                    println("                   Effective Precision: ", preceff)
                end
            end
        end
    end
end