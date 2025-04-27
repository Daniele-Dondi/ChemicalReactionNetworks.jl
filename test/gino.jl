using ChemicalReactionNetworks

# example taken from Sarkar paper
reactions = [Reaction([1,1],[1,1],[2],[1],1.0,1.0),
            Reaction([1,2],[1,1],[3],[1],1.0,1.0),
            Reaction([1,3],[1,1],[4],[1],1.0,1.0),
            Reaction([1,3],[1,1],[2],[2],2.0,1.0),
            Reaction([1,4],[1,1],[2,3],[1,1],1.0,1.0),
            Reaction([2,4],[1,1],[3],[2],1.0,1.0)]
chem_potentials=ones(4)
fe=free_energies(chem_potentials, reactions)
kf,kr = reaction_rates(fe, ones(6), 1.0)
for i=1:6
    reactions[i].kf = kf[i]
    reactions[i].kr = kr[i]
end 
z=equilibrium_state(ReactionNetwork(reactions))
print(z)
dz=mass_action(reactions, z)
print(dz)
@assert all(dz .< 1e-9)   # the network should admit an equilibrium state

