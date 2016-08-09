#!/usr/bin/perl

my $input; # input file name

use strict;
use warnings;
use Data::Dumper;
use Switch;
use List::Util 1.33 'any';
use List::Util 'first';
use Math::BigFloat;

# buffer size (in number of events)
my $buffer_size = 10000;

# use dictionary to keep track of states and events
my %states;
my %events;

# use dictionary to keep track of translated and ignored events
my %translated_events;
my %ignored_events;

my $number_of_tasks;

# options
my $prv_time_unit = "ns";

# store the app end time
my $app_end_time = "0";

# store the pjdump container
my @pjdump_containers;

# store the pjdump link
my @pjdump_links;

# store the pjdump states
my @pjdump_states;

# store the pjdump events
my @pjdump_events;

# store the current state value for each task
my @task_current_state_value;


# store communicator id and size
my %communicators;

my @task_states_buffer;
my @task_events_buffer;
my @task_comms_buffer;

my $power_reference = 286.087E-3; # in flop/mus



sub getTime {
    my($timeValue) = @_;

    if ($prv_time_unit eq 'ns') {
	return ($timeValue / 1000000000);
    }

    return $timeValue;
}

sub pjdump {

    foreach (@pjdump_containers) {
        my $action = $_;
	my $parent = $action->{"parent"};
	my $type = $action->{"type"};
        my $startTime = getTime($action->{"startTime"});
        my $endTime = getTime($app_end_time);
        my $duration = getTime($app_end_time - $action->{"startTime"});
        my $name = $action->{"name"};
	print "Container, $parent, $type, $startTime, $endTime, $duration, $name\n";
    }

    foreach (@pjdump_links) {
        my $action = $_;
	my $container = $action->{"container"};
	my $type = $action->{"type"};
        my $startTime = getTime($action->{"startTime"});
        my $endTime = getTime($action->{"endTime"});
        my $duration = getTime($action->{"duration"});
        my $value = $action->{"value"};
        my $startContainer = $action->{"startContainer"};
        my $endContainer = $action->{"endContainer"};
	print "Link, $container, $type, $startTime, $endTime, $duration, $value, $startContainer, $endContainer\n";
    }

    foreach (@pjdump_states) {
        my $state = $_;
	my $container = $state->{"container"};
	my $type = $state->{"type"};
        my $startTime = getTime($state->{"startTime"});
        my $endTime = getTime($state->{"endTime"});
        my $duration = getTime($state->{"duration"});
        my $imbrication = $state->{"imbrication"};
        my $value = $state->{"value"};
	print "State, $container, $type, $startTime, $endTime, $duration, $imbrication, $value\n";
    }

    foreach (@pjdump_events) {
        my $action = $_;
	my $container = $action->{"container"};
	my $name = $action->{"name"};
        my $time = getTime($action->{"time"});
        my $value = $action->{"value"};
	print "Event, $container, $name, $time, $value\n";
    }

    # clear data structures after pjdump
    @pjdump_containers = ();
    @pjdump_links = ();
    @pjdump_states = ();
    @pjdump_events = ();
}

sub main {
    my($arg);
    
    while(defined($arg = shift(@ARGV))) {
        for ($arg) {
            if (/^-i$/) { $input = shift(@ARGV); last; }
            print "unrecognized argument '$arg'\n";
        }
    }
    if(!defined($input) || $input eq "") { die "No valid input file provided.\n"; }
    
    parse_pcf($input.".pcf");
    parse_prv_lucas($input.".prv");

    # pjdump();
}

my %mpi_call_parameters = (
    "send size" => "50100001",
    "recv size" => "50100002",
    "root" => "50100003",
    "communicator" => "50100004",
    );

my @mpi_calls = (
    "MPI_Finalize",       #
    "MPI_Init",           #
    "MPI_Send",           #
    "MPI_Recv",           #
    "MPI_Isend",          #
    "MPI_Irecv",          #
    "MPI_Wait",           #
    "MPI_Waitall",        #
    "MPI_Bcast",          #
    "MPI_Reduce",         #
    "MPI_Allreduce",      #
    "MPI_Barrier",        #
    "MPI_Gather",         #
    "MPI_Allgather",      #
    "MPI_Alltoall",       #
    "MPI_Gatherv",        #
    "MPI_Allgatherv",     #
    "MPI_Reduce_scatter", #
    "MPI_Alltoallv"
    );


# search for a MPI call in the event's parameters
# in all the cases I have seen, the event type and value are the first
# numbers in the event's parameter list, however we are not making this
# assumption. Instead, we look at all parameters and search for the one
# that is encoding the MPI call
sub extract_mpi_call {
    my %event_info = @_;
    
    # search for a MPI call in the event's parameters
    foreach my $key (keys %event_info) {
	if(defined($events{$key})) {
	    if(defined($events{$key}{value}{$event_info{$key}})) {
		my $event_name = $events{$key}{value}{$event_info{$key}};

                if (any { /^$event_name$/ } @mpi_calls){
		    $translated_events{$event_name} = 1;
		    return $event_name;
		}
		else {
		    $ignored_events{$event_name} = 1;
		}
	    }
	}
    }
    return "None";
}


sub parse_prv_lucas {
    my($prv) = @_; # get arguments
    open(INPUT, $prv) or die "Cannot open $prv. $!";

    # check if header is valid, we should get something like #Paraver (dd/mm/yy at hh:m):ftime:0:nAppl:applicationList[:applicationList]
    my $line = <INPUT>;
    chomp $line;
    $line =~ /^\#Paraver / or die "Invalid header '$line'\n";
    my $header = $line;
    $header =~ s/^[^:\(]*\([^\)]*\):// or die "Invalid header '$line'\n";
    $header =~ s/(\d+):(\d+)([^\(\d])/$1\_$2$3/g;
    $header =~ s/,\d+$//g;
    my($max_duration, $resource, $number_of_apps, @app_info_list) = split(/:/, $header);
    $max_duration =~ s/_.*$//g;
    $resource =~ /^(.*)\((.*)\)$/ or die "Invalid resource description '$resource'\n";
    my($number_of_nodes, $node_cpu_count) = ($1, $2);
    $number_of_apps == 1 or die "Only one application can be handled at the moment\n";
    my @node_cpu_count = split(/,/, $node_cpu_count);

    # parse app info
    foreach my $app (1..$number_of_apps) {
        $app_info_list[$app - 1] =~ /^(.*)\((.*)\)$/ or die "Invalid application description\n";
	my $task_info;
        ($number_of_tasks, $task_info) = ($1, $2);
        my(@task_info_list) = split(/,/, $task_info);

        # create app container
        my %container;
        $container{"parent"} = "0";
        $container{"type"} = "APP";
        $container{"name"} = "0";
        $container{"startTime"} = "0";
        push @pjdump_containers, \%container;

        # initialize current state value for each task
	@task_current_state_value = (0) x $number_of_tasks;;

        # initiate an empty event buffer for each task
	@task_events_buffer = (0);
        foreach my $task (1..$number_of_tasks) {
            my($number_of_threads, $node) = split(/_/, $task_info_list[$task - 1]);

            # create task container
            my %container;
            $container{"parent"} = "0";
            $container{"type"} = "TASK";
            $container{"name"} = ($task - 1);
            $container{"startTime"} = "0";
            push @pjdump_containers, \%container;
        }
    }

    # LUCAS version LUCAS
    # start reading records
    my $event_counter = 0;
    while(defined($line=<INPUT>)) {
        chomp $line;

        # state records are in the format 1:cpu:appl:task:thread:begin_time:end_time:state_id
        if($line =~ /^1/) {
            my($record, $cpu, $appli, $task, $thread, $begin_time, $end_time, $state_id) = split(/:/, $line);
	    my $state_name = $states{$state_id}{name};

	    # convert the time values to high precision numbers
	    #$begin_time = Math::BigFloat->new($begin_time);
	    #$end_time = Math::BigFloat->new($end_time);

            # # create state
            # my %state;
            # $state{"container"} = "rank-" . ($task - 1);
            # $state{"type"} = uc($state_name);
            # $state{"startTime"} = $begin_time;
	    # $state{"endTime"} = $end_time;
            # $state{"duration"} = $end_time - $begin_time;
            # $state{"imbrication"} = "0";
            # $state{"value"} = uc($state_name);

	    print "rank-" . ($task - 1) . ", $begin_time, $end_time, \"$state_name\"\n";
	    
            #push @pjdump_states, \%state;

            # update the reference to the value of the current state of this task
            #$task_current_state_value[$task - 1] = \$state{"value"};

	    # if ($end_time + 0 > $app_end_time + 0) {
	    # 	$app_end_time = $end_time;
            # }
        }
        
	# # event records are in the format 2:cpu:appl:task:thread:time:event_type:event_value
        # elsif($line =~ /^2/){
	#     continue;
	#     my($record, $cpu, $appli, $task, $thread, $time, %event_list) = split(/:/, $line);
	#     my $mpi_call = extract_mpi_call(%event_list);

	#     # convert the time value to a high precision number
	#     $time = Math::BigFloat->new($time);

	#     # if event is a MPI call, get the MPI call parameters and add_event_entry($task, %parameters);
	#     if($mpi_call ne "None") {

        #         # create event
        #         my %event;
        #         $event{"container"} = ($task - 1);
        #         $event{"name"} = "MPI_CALL";
        #         $event{"time"} = $time;
        #         $event{"value"} = uc($mpi_call);
        #         push @pjdump_events, \%event;

        #         #my $a = $task_current_state_value[$task - 1];
        #         #$$a = uc($mpi_call);
	#     }     
        # }

        # # communication records are in the format 3:cpu_send:ptask_send:task_send:thread_send:logical_time_send:actual_time_send:cpu_recv:ptask_recv:task_recv:thread_recv:logical_time_recv:actual_time_recv:size:tag
	# elsif($line =~ /^3/) {
	#     continue;
        #    my($record, $cpu_send, $ptask_send, $task_send, $thread_send, $ltime_send, $atime_send, $cpu_recv, $ptask_recv, $task_recv, $thread_recv, $ltime_recv, $atime_recv, $size, $tag) = split(/:/, $line);

	#    # convert the time values to high precision numbers
	#    $atime_send = Math::BigFloat->new($atime_send);
	#    $atime_recv = Math::BigFloat->new($atime_recv);

        #    # create link
        #    my %link;
        #    $link{"container"} = "0";
        #    $link{"type"} = "LINK";
        #    $link{"startTime"} = $atime_send;
	#    $link{"endTime"} = $atime_recv;
        #    $link{"duration"} = $atime_recv - $atime_send;
        #    $link{"value"} = "LINK";
        #    $link{"startContainer"} = ($task_send - 1);
        #    $link{"endContainer"} = ($task_recv - 1);
        #    push @pjdump_links, \%link;
        # }

	# $event_counter++;
	# if ($event_counter > $buffer_size) {
	#     pjdump();
	#     $event_counter = 0;
	#     exit();
	# }
    }

    return;
}



sub parse_pcf {
    my($pcf) = shift; # get first argument
    my $line;

    open(INPUT, $pcf) or die "Cannot open $pcf. $!";
    while(defined($line=<INPUT>)) {
        chomp $line; # remove new line
        if($line =~ /^STATES$/) {
            while((defined($line=<INPUT>)) && ($line =~ /^(\d+)\s+(.*)/g)) {
                $states{$1}{name} = $2;
		$states{$1}{used} = 0;
            }
        }

        if($line =~ /^EVENT_TYPE$/) {
	    my $id;
            while($line=<INPUT>) { # read event
                if($line =~ /VALUES/g) {
		    while((defined($line=<INPUT>)) && ($line =~ /^(\d+)\s+(.*)/g)) { # read event values
			$events{$id}{value}{$1} = $2;
		    }
		    last;
		}
                $line =~ /[\d]\s+(\d+)\s+(.*)/g or next;
                $id = $1;
                $events{$id}{type} = $2;
		$events{$id}{used} = 0;
            }
        }
    }
}


main();
