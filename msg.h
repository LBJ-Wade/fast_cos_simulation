#ifndef MSG_H
#define MSG_H 1

enum LogLevel {verbose, debug, normal, info, warn, error, fatal, silent};

//void msg_init(void);
void msg_set_loglevel(const enum LogLevel level);
void msg_printf(const enum LogLevel level, const char *fmt, ...);
void msg_abort(const char *fmt, ...);

#endif
