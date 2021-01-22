#include <stdio.h>
#ifdef HAVE_PROJ_H
#include <proj.h>
int main() {
    printf("%d%d%d\n", PROJ_VERSION_MAJOR, PROJ_VERSION_MINOR, PROJ_VERSION_PATCH);
    return 0;
}
#else
#include <proj_api.h>
int main() {
    printf("%d\n", PJ_VERSION);
    return 0;
}
#endif
