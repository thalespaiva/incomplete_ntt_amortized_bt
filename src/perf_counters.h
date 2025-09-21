// THIS FILE IS BASED ON:
// https://learn.arm.com/learning-paths/servers-and-cloud-computing/arm_pmu/perf_event_open/

#pragma once

#include <linux/perf_event.h> /* Definition of PERF_* constants */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/ioctl.h>
#include <sys/syscall.h> /* Definition of SYS_* constants */
#include <unistd.h>
#include <inttypes.h>
#define PERF_COUNTERS_TOTAL_EVENTS 2

// The function to counting through (called in main)
static void code_to_measure(){
  int sum = 0;
  for(int i = 0; i < 1000000000; ++i){
    sum += 1;
  }
}

// Executes perf_event_open syscall and makes sure it is successful or exit
static long perf_event_open(struct perf_event_attr *hw_event, pid_t pid, int cpu, int group_fd, unsigned long flags){
  int fd;
  fd = syscall(SYS_perf_event_open, hw_event, pid, cpu, group_fd, flags);
  if (fd == -1) {
    fprintf(stderr, "Error creating event");
    exit(EXIT_FAILURE);
  }

  return fd;
}

// Helper function to setup a perf event structure (perf_event_attr; see man perf_open_event)
static void configure_event(struct perf_event_attr *pe, uint32_t type, uint64_t config){
  memset(pe, 0, sizeof(struct perf_event_attr));
  pe->type = type;
  pe->size = sizeof(struct perf_event_attr);
  pe->config = config;
  pe->read_format = PERF_FORMAT_GROUP | PERF_FORMAT_ID;
  pe->disabled = 1;
  pe->exclude_kernel = 1;
  pe->exclude_hv = 1;
}

// Format of event data to read
// Note: This format changes depending on perf_event_attr.read_format
// See `man perf_event_open` to understand how this structure can be different depending on event config
// This read_format structure corresponds to when PERF_FORMAT_GROUP & PERF_FORMAT_ID are set
struct read_format {
  uint64_t nr;
  struct {
    uint64_t value;
    uint64_t id;
  } values[PERF_COUNTERS_TOTAL_EVENTS];
};

int fd[PERF_COUNTERS_TOTAL_EVENTS];  // fd[0] will be the group leader file descriptor
size_t id[PERF_COUNTERS_TOTAL_EVENTS];  // event ids for file descriptors
uint64_t pe_val[PERF_COUNTERS_TOTAL_EVENTS]; // Counter value array corresponding to fd/id array.
struct perf_event_attr pe[PERF_COUNTERS_TOTAL_EVENTS];  // Configuration structure for perf events (see man perf_event_open)
struct read_format counter_results;


static void perf_counters_start() {
  configure_event(&pe[0], PERF_TYPE_HARDWARE, PERF_COUNT_HW_CACHE_MISSES);
  configure_event(&pe[1], PERF_TYPE_HARDWARE, PERF_COUNT_HW_CACHE_REFERENCES);
  // configure_event(&pe[3], PERF_TYPE_HARDWARE, PERF_COUNT_HW_STALLED_CYCLES_BACKEND);
  // configure_event(&pe[4], PERF_TYPE_RAW, 0x70);  // Count of speculative loads (see Arm PMU docs)
  // configure_event(&pe[5], PERF_TYPE_RAW, 0x71);  // Count of speculative stores (see Arm PMU docs)

  // Create event group leader
  fd[0] = perf_event_open(&pe[0], 0, -1, -1, 0);
  ioctl(fd[0], PERF_EVENT_IOC_ID, &id[0]);
  // Let's create the rest of the events while using fd[0] as the group leader
  for(int i = 1; i < PERF_COUNTERS_TOTAL_EVENTS; i++){
    fd[i] = perf_event_open(&pe[i], 0, -1, fd[0], 0);
    ioctl(fd[i], PERF_EVENT_IOC_ID, &id[i]);
  }

  // Reset counters and start counting; Since fd[0] is leader, this resets and enables all counters
  // PERF_IOC_FLAG_GROUP required for the ioctl to act on the group of file descriptors
  ioctl(fd[0], PERF_EVENT_IOC_RESET, PERF_IOC_FLAG_GROUP);
  ioctl(fd[0], PERF_EVENT_IOC_ENABLE, PERF_IOC_FLAG_GROUP);
}


static void perf_counters_stop() {
  // Stop all counters
  ioctl(fd[0], PERF_EVENT_IOC_DISABLE, PERF_IOC_FLAG_GROUP);

  // Read the group of counters and print result
  read(fd[0], &counter_results, sizeof(struct read_format));
  // printf("Num events captured: %" PRIu64 "\n", counter_results.nr);
  for(size_t i = 0; i < counter_results.nr; i++) {
    for(size_t j = 0; j < PERF_COUNTERS_TOTAL_EVENTS ;j++){
      if(counter_results.values[i].id == id[j]){
        pe_val[i] = counter_results.values[i].value;
      }
    }
  }
  // printf("PERF_COUNT_HW_CACHE_MISSES: %" PRIu64 "\n", pe_val[0]);
  // printf("PERF_COUNT_HW_CACHE_REFERENCES: %" PRIu64 "\n", pe_val[1]);

  printf("cache miss rate: %ld / %ld = %.4lf\n", pe_val[0], pe_val[1], (double) pe_val[0] / pe_val[1]);

  // Close counter file descriptors
  for(int i = 0; i < PERF_COUNTERS_TOTAL_EVENTS; i++){
    close(fd[i]);
  }
}