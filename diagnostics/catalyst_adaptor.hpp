#ifndef CATALYST_ADAPTOR_H_
#define CATALYST_ADAPTOR_H_

#if USE_CATALYST

extern "C" {
void catalyst_adaptor();

void catalyst_adaptor_initialise(char *a_catalyst_script);

void catalyst_adaptor_finalise();
}

#endif /* USE_CATALYST */
#endif