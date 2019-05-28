function aE_generateNWB(dirs,xlsdata)
savefield=false;
if isfield(dirs,'NWBdir')
    %% calculating cell coordinates 
    cellcoord_lat=[xlsdata.Cranio_Lat]+(([xlsdata.locationX]-[xlsdata.Cranio_center_X])/1000).*[xlsdata.Lateral_dir_X]; % works only when animal's lateral axis is sideways to the stage
    cellcoord_AP=[xlsdata.Cranio_AP]+(([xlsdata.locationy]-[xlsdata.Cranio_center_Y])/1000).*[xlsdata.Rostral_dir_Y]; % works only when the animal's a-p axis is parallel to the stage
    for i=1:length(xlsdata)
        xlsdata(i).cellcoord_lat=cellcoord_lat(i);
        xlsdata(i).cellcoord_AP=cellcoord_AP(i);
    end
    %%
    for xlsidx=length(xlsdata):-1:1
        if xlsdata(xlsidx).field==0 & xlsdata(xlsidx).juxta==0 % runs only for intracellular recordings
            ID=xlsdata(xlsidx).ID;
            load([dirs.bridgeddir,ID],'lightdata','bridgeddata','stimdata');
            load([dirs.rawexporteddir,ID],'rawdata');
            source_file = [mfilename() '.m'];
            [~, source_script, ~] = fileparts(source_file);
            nwb = nwbfile();
            nwb.identifier = ID;
            nwb.general_source_script = source_script;
            nwb.general_source_script_file_name = source_file;
            nwb.general_lab = 'Gabor Tamas - MTA-SZTE Research Group for Cortical Microcircuits';
%             nwb.general_keywords = {''};
            nwb.general_institution = ['University of Szeged,'...
                ' Department of Physiology, Anatomy and Neuroscience, Szeged, Hungary'];
            if any(strfind(ID,'og'))
                experimenter='Gáspár Oláh';
            elseif any(strfind(ID,'rm'))
                experimenter='Márton Rózsa';
            elseif any(strfind(ID,'tm'))
                experimenter='Martin Tóth';
            elseif any(strfind(ID,'mn'))
                experimenter='Norbert Mihut';
            elseif any(strfind(ID,'kb'))
                experimenter='Balázs Kovács';
            elseif any(strfind(ID,'mg'))
                experimenter='Gábor Molnár';
            end
            nwb.general_experimenter=experimenter;
            %             nwb.general_related_publications = [''];
            %             nwb.general_stimulus = 'sound';
            %             nwb.general_protocol = '';
            %             nwb.general_surgery = [''];
            nwb.session_description = [ID,' - in vivo recording'];
            nwb.session_start_time = datetime([ID(1:6),datestr(seconds(lightdata(1).realtime),'HHMMSSFFF')],'InputFormat', 'yyMMddHHmmssSSS');
            nwb.timestamps_reference_time = nwb.session_start_time;
            nwb.general_subject = types.core.Subject(...
                'species', xlsdata(xlsidx).species, ...
                'age', [num2str(xlsdata(xlsidx).age),' days postnatal'], ...
                'description', ['anesthesia: ', xlsdata(xlsidx).anaesthesia]);
            %% adding whole cell electrode
            ic_device_name = rawdata(1).AmplifierID{(rawdata(1).tracenumber)};
            ic_elec_name = 'Whole cell patch clamp electrode';
            nwb.general_devices.set(ic_device_name, types.core.Device());
            ic_device_link = types.untyped.SoftLink(['/general/devices/' ic_device_name]);
            
            ic_elec = types.core.IntracellularElectrode( ...
                'device', ic_device_link, ...
                'description', 'Whole cell patch clamp electrode pulled with DMZ universal puller, 4-6 MOhm');
            
            nwb.general_intracellular_ephys.set(ic_elec_name, ic_elec);
            ic_elec_link = types.untyped.SoftLink(['/general/intracellular_ephys/' ic_elec_name]);
            
            %% checking for field and adding electrode
            if isfield(xlsdata,'field') & savefield
                fieldxlsnum=find(strcmp(xlsdata(xlsidx).HEKAfname,{xlsdata.HEKAfname}) & [xlsdata.field]==1);
            else
                fieldxlsnum=[];
            end
            if length(fieldxlsnum)>1
                disp('error, more than 1 field files found.. user should select..')
                ok=0;
                Selection=1;
                while ok==0 | length(Selection)~=1
                    [Selection,ok] = listdlg('PromptString','Select the corresponding field file:','ListString',{xlsdata(fieldxlsnum).ID});
                end
                fieldxlsnum=fieldxlsnum(Selection);
            end
            if ~isempty(fieldxlsnum)
                field_Z=xlsdata(fieldxlsnum).locationz;
                if isnan(xlsdata(fieldxlsnum).cellcoord_AP)
                    field_X=xlsdata(fieldxlsnum).locationX;
                    field_Y=xlsdata(fieldxlsnum).locationy;
                    location_description='relative to the recording rig';
                else
                    field_X=xlsdata(fieldxlsnum).cellcoord_lat;
                    field_Y=xlsdata(fieldxlsnum).cellcoord_AP;
                    location_description='relative to bregma (x = lateral, y = anterior-posterior)';
                end
                
                
                ID_field=xlsdata(fieldxlsnum).ID;
                rawdata_field=load([dirs.rawexporteddir,ID_field],'rawdata');
                rawdata_field=rawdata_field.rawdata;
                
                field_device_name = rawdata_field(1).AmplifierID{(rawdata_field(1).tracenumber)};
                %                 field_elec_name = ['Extracellular glass electrode'];
                nwb.general_devices.set(field_device_name, types.core.Device());
                field_device_link = types.untyped.SoftLink(['/general/devices/' field_device_name]);
                egroup = types.core.ElectrodeGroup(...
                    'description', 'Extracellular glass electrode',...
                    'location', xlsdata(xlsidx).corticalregion,...
                    'device', field_device_link);
                nwb.general_extracellular_ephys.set(field_device_name, egroup);
                
                %% this doesn't work for me...
                dtColNames = {'x', 'y', 'z'};
                dynTable = types.core.DynamicTable(...
                    'colnames', dtColNames,...
                    'description', 'Electrodes',...
                    'id', types.core.ElementIdentifiers('data', 1),...
                    'x', types.core.VectorData('data', field_X, 'description', ['the x coordinate of the electrode ', location_description]),...
                    'y', types.core.VectorData('data', field_Y, 'description', ['the y coordinate of the electrode ', location_description]),...
                    'z', types.core.VectorData('data', field_Z, 'description','the z coordinate of the channel relative to the surface of the brain'));
                nwb.general_extracellular_ephys.set('electrodes', dynTable);
            end
            %% adding every sweep
            for sweep = 1:length(bridgeddata)
                %% adding voltage response
                data = single(bridgeddata(sweep).y); % using single precision for compression
                response_name = ['Whole cell - Sweep ',num2str(sweep)];
                nwb.acquisition.set(response_name, ...
                    types.core.CurrentClampSeries( ...
                    'bias_current', [], ... % Unit: Amp
                    'bridge_balance', stimdata(sweep).RS, ... % Unit: Ohm
                    'capacitance_compensation', [], ... % Unit: Farad
                    'data', data, ...
                    'data_unit', 'V', ...
                    'electrode', ic_elec_link, ...
                    'stimulus_description', 'Intracellular whole cell recording', ...
                    'sweep_number', uint64(sweep), ...
                    'starting_time', bridgeddata(sweep).realtime, ...
                    'starting_time_rate', round(1/bridgeddata(sweep).si)));
                
                %% adding stimulus
                data = single(stimdata(sweep).y); % using single precision for compression
                stimulus_name = ['Whole cell - Sweep ',num2str(sweep)];
                nwb.stimulus_presentation.set(stimulus_name, ...
                    types.core.CurrentClampStimulusSeries( ...
                    'data', data, ...
                    'data_unit', 'A', ...
                    'electrode', ic_elec_link, ...
                    'stimulus_description', 'Intracellular whole cell recording', ...
                    'sweep_number', uint64(sweep), ...
                    'starting_time', bridgeddata(sweep).realtime, ...
                    'starting_time_rate', round(1/bridgeddata(sweep).si)));
                
                %% adding extracellular recording
                if ~isempty(fieldxlsnum)
                    fieldsweepnum=find([rawdata_field.realtime]==bridgeddata(sweep).realtime);
                    if ~isempty(fieldsweepnum)
            
                        tablereg = types.core.DynamicTableRegion(...
                            'description','Electrode for this Electrical Series',...
                            'table',types.untyped.ObjectView('/general/extracellular_ephys/electrodes'),...
                            'data',[1 1]);
                        %%
                        data = single(rawdata_field(fieldsweepnum).y); % using single precision for compression
                        field_name = ['Raw extracellular recording - sweep ',num2str(sweep)];
                        nwb.acquisition.set(field_name, ...
                            types.core.ElectricalSeries( ...
                            'data', data, ...
                            'data_unit', 'V', ...% %
                            'electrodes', tablereg,... this was skipped..
                            'starting_time', bridgeddata(sweep).realtime, ...
                            'starting_time_rate', round(1/bridgeddata(sweep).si)));
                        
                    end
                end
                
            end

            nwbExport(nwb, [dirs.NWBdir,ID,'.nwb']);
            disp([ID,' exported to NWB'])
%             return
        end
    end
end