package com.wyrm.jscheduler.Entity;

import jakarta.persistence.Column;
import jakarta.persistence.Entity;
import lombok.AllArgsConstructor;
import lombok.Data;
import lombok.NoArgsConstructor;
import jakarta.persistence.Id;
import jakarta.persistence.Table;
import org.json.JSONObject;
import org.hibernate.annotations.Type;
import org.hibernate.annotations.JdbcTypeCode;
import org.hibernate.type.SqlTypes;
import java.time.Instant;
import java.util.UUID;
import io.hypersistence.utils.hibernate.type.json.JsonBinaryType;

@NoArgsConstructor
@AllArgsConstructor
@Data
@Entity
@Table(name="task")
public  class Task
{
    @Id
    private UUID iD;
    private String type;
    private JobStatus status;
    @Type(JsonBinaryType.class)
    @JdbcTypeCode(SqlTypes.JSON)
    @Column(name="data", columnDefinition = "jsonb")
    private JSONObject data;
    private Instant timestamp;
    private int attempts;
    public enum JobStatus
    {
        PENDING,
        RETRYING,
        FAILED,
        COMPLETED
    }
}