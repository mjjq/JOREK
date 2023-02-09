#include "catalyst_adaptor.hpp"

#if USE_CATALYST
#include "catalyst.hpp"
#include <iostream>

extern "C" {
void catalyst_adaptor() {}

void catalyst_adaptor_initialise(char *a_catalyst_script) {
  std::cout << "The Catalyst script passed to me is:\n"
            << a_catalyst_script << std::endl;
  conduit_cpp::Node node;

  // Pass script to Catalyst
  node["catalyst/scripts/script0"].set_string(a_catalyst_script);

  // Initialize Catalyst
  catalyst_status err = catalyst_initialize(conduit_cpp::c_node(&node));
  if (err != catalyst_status_ok) {
    std::cerr << "Failed to initialize Catalyst: " << err << std::endl;
  }
}

void catalyst_adaptor_finalise() {
  conduit_cpp::Node node;

  // Finalize Catalyst
  catalyst_status err = catalyst_finalize(conduit_cpp::c_node(&node));
  if (err != catalyst_status_ok) {
    std::cerr << "Failed to finalize Catalyst : " << err << std::endl;
  }
}
}

#endif