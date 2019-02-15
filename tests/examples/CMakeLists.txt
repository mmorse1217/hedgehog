
set(Hedgehog_test_examples_src
    matvec_examples.cpp solver_examples.cpp
)

PREPEND(Hedgehog_test_examples_src)

set(Hedgehog_test_examples ${Hedgehog_test_examples_src} PARENT_SCOPE)
