#ifndef FUSION_PERFCOST_H
#define FUSION_PERFCOST_H

#include "perfmon_wrapper.h"
#include "cost_model.h"

extern struct perf_data *global_perf_data;
extern bool global_perf_measuring;

#define PERF_FUNC_ENTER {COST_FUNC_ENTER; if(global_perf_measuring){perf_start(global_perf_data);} }
#define PERF_FUNC_EXIT { COST_FUNC_EXIT; if(global_perf_measuring){perf_stop(global_perf_data);} }
#define PERF_FUNC_RESET {perf_reset(global_perf_data); COST_MODEL_RESET;}

#endif // FUSION_PERFCOST_H
