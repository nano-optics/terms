The authors describe in this manuscript (or rather, its supplementary information), the application of the Born-Kuhn model of coupled oscillators to the description of chiroptical effects between an achiral plasmonic nanoparticle and a chiral nanocrystal. 

The main manuscript "sells" the model by showing simulations of several configurations, tuning the different parameters of the model with arbitrary values to illustrate the effect of weak or strong coupling between resonators, and the effect of detuning. Fully numerical simulations using Finite elements are then used in separate plots to show that these quantitative behaviours can indeed be reproduced with fully explicit numerical simulations.

I believe the model has some value, as it provides physical insights and is computationally efficient. The results seem plausible to me, qualitatively, and the physics is probably sound (the Born-Kuhn model is certainly well-accepted). However, I have strong reservations regarding the presentation, structure, as well as the novelty and breadth of this study. I do not think that in its present form it is suitable for ACS Nano, but perhaps with some changes discussed below it could be of interest in a more technical journal of the ACS family.

## Presentation and structure

- The manuscript is constantly referring to the S.I., and it's only when I read the S.I. that I knew what simulations were actually done, and what the model actually is.

- A lot of the manuscript is trying to "sell" the model, insisting on the quality of its results etc. I think this is not helpful for a scientific publication; I would expect a more critical analysis of the model, its pros as well as its limitations, a thorough discussion of its assumptions and the approximations used, a thorough comparison with, and discussion of, "competing" methods such as the coupled-dipole method, multipolar models, Lagrangian models, effective medium models, etc. I saw none of this in the manuscript. And I was not "sold" on this model.

- The FEM results are presented as confirmation of the model -- but in separate figures, with absolutely no attempt at a quantitative comparison. 

- The model is not described very well: what is even calculated? Absorption? Extinction? Why is everything in arbitrary units? Why are all the parameters arbitrary? Could they be obtained for a specific structure? How? If one has to use arbitrary coupling parameters, then how is this model useful (beside providing a conceptual insight into the physics, which I agree is useful -- as was the Born-Kuhn model all along).


- While the introduction was relatively clear, I found the rest of the manuscript quickly becoming much harder to follow (unclear), and not just because most of the explanations are referring to the S.I.


## Breadth


- What is the incident field? Does it come along a specific direction? CD is often orientation-averaged, can the model do that? How do the responses vary with angle?

- The plasmonic particle ("oscillator" sic) is said to flip the phase of the local electric field (which seems reasonable), but surely there's more to it: what about the vector aspect? The field near a metal surface is very different from the field in a homogeneous dielectric. Again the complete lack of information about the incident field is quite problematic, and the effect of a nanoparticle on its near-field distribution should be discussed in depth (including, I suspect, many limitations of the model discussed here).

- I'm not sure much of this model is very new. Giessen and co-workers described the "plasmonic Born-Kuhn" model many years ago. Govorov and co-workers described the coupling between a plasmonic nanoparticle and a chiral unit over a decade ago. 

- From the S.I. I finally gathered that retardation is neglected. This seems like a very important assumption and limitation. Note that it is not in fact consistent, as the incident field's retardation is kept. (Removing it would make the CD vanish...) Neglecting the retardation between oscillators will distort the picture, in general, except for very small structures (again this has been known for at least a decade). I would expect at least to see a detailed study of this approximation, its range of validity, the limitations, etc.

- I don't find the "zero-order" treatment of numerical error in FEM simulations very convincing. It is plausible, but very context dependent. If a structure has a complicated mesh, if the structure breaks obvious symmetries, etc. all kinds of interactions could happen where the interplay of artifical chirality of the mesh and of the idealised system could easily mix and make things more complicated.

