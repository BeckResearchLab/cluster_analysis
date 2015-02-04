SET @instance_pid = NULL;
SET @instance_time = NULL;
SET @session_id = NULL;
SET @clusterLike = NULL;
SET @myClusterLike = NULL;

SELECT * FROM (

SELECT * ,
		IF(`instance.pid` = @instance_pid AND
			`instance.time` = @instance_time AND
			`session.id` = @session_id,
				IF(clusterLike != @clusterLike, "clusterLike",
					IF(myClusterLike != @myClusterLike, "myClusterLike",
						""
					)
				) AS liked_how,
		@instance_pid := `instance.pid`,
		@instance_time := `instance.time`,
		@session_id := `session.id`,
		@clusterLike := clusterLike
		@myClusterLike := myClusterLike
	FROM log
	ORDER BY `instance.time`, `instance.pid`, `session.id`, `session.counter`

) AS ordered_log 
	WHERE NOT(ISNULL(ordered_log.liked_how))
;
