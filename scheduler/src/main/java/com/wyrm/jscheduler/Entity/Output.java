package com.wyrm.jscheduler.Entity;

import io.hypersistence.utils.hibernate.type.json.JsonBinaryType;
import jakarta.persistence.Column;
import jakarta.persistence.Entity;
import lombok.AllArgsConstructor;
import lombok.Data;
import lombok.NoArgsConstructor;
import jakarta.persistence.Id;
import jakarta.persistence.Table;
import org.hibernate.annotations.JdbcTypeCode;
import org.hibernate.annotations.Type;
import org.hibernate.type.SqlTypes;
import org.json.JSONObject;

import java.time.Instant;
import java.util.UUID;

@NoArgsConstructor
@AllArgsConstructor
@Data
@Entity
@Table(name="output")
public  class Output
{
    @Id
    private UUID ID;
    private String type;
    @Type(JsonBinaryType.class)
    @JdbcTypeCode(SqlTypes.JSON)
    @Column(name="data", columnDefinition = "json")
    private JSONObject results;
    private int attempts;
    public enum JobStatus
    {
        PENDING,
        RETRYING,
        FAILED,
        COMPLETED
    }
}