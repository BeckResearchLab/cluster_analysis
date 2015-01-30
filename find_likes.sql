SET @instance_pid = NULL;
SET @instance_time = NULL;
SET @session_id = NULL;
SET @myClusterLike = NULL;

SELECT * FROM (

SELECT * ,
		IF(`instance.pid` = @instance_pid AND
			`instance.time` = @instance_time AND
			`session.id` = @session_id AND 
			myClusterLike != @myClusterLike, 1, 0) AS liked,
		@instance_pid := `instance.pid`,
		@instance_time := `instance.time`,
		@session_id := `session.id`,
		@myClusterLike := myClusterLike
	FROM log
	ORDER BY `instance.time`, `instance.pid`, `session.id`, `session.counter`

) AS ordered_log 
	WHERE ordered_log.liked = 1
;
