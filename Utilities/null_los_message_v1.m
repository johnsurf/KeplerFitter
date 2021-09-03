function [message] = null_NER_message()

% % Header Row Data Fields
message.trackID   =[];
message.blockID   = [];
message.availTime = [];

% Message Block, row 2
message.timeUTC   = [];
message.pad1      = [];
message.satNum    = [];
message.losVec    = zeros(1,3);
message.satVec    = zeros(1,3);
message.radioInt  = [];
message.band      = [];
message.pad2      = [];

end