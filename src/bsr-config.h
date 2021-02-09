#ifndef BSR_CONFIG_H
#define BSR_CONFIG_H

void initConfig(bsr_config_t *bsr_config);
int loadConfigFromFile(bsr_config_t *bsr_config);
int loadConfigFromQueryString(bsr_config_t *bsr_config, char *query_string);
int validateCGIOptions(bsr_config_t *bsr_config);

#endif // BSR_CONFIG_H
