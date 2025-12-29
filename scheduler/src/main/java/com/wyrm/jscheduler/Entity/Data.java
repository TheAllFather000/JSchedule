package com.wyrm.jscheduler.Entity;
import java.util.UUID;

public record Data(
        UUID id,
        String output,
    double processing_time
){};