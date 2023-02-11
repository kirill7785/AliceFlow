unit Unitdefineindentcabinet;
  // Формирование размеров кабинета по существующему наполнению
  // расчётной модели блоками.
interface

uses
  Winapi.Windows, Winapi.Messages, System.SysUtils, System.Variants, System.Classes, Vcl.Graphics,
  Vcl.Controls, Vcl.Forms, Vcl.Dialogs, Vcl.StdCtrls, Vcl.ExtCtrls,
  Vcl.Imaging.pngimage;

type
  TFormcabinetindent = class(TForm)
    Image1: TImage;
    RadioGroup1: TRadioGroup;
    GroupBox1: TGroupBox;
    ButtonNext: TButton;
    PanelH_only: TPanel;
    Label1: TLabel;
    Edithonly: TEdit;
    Labelunion: TLabel;
    Panelhxyz: TPanel;
    Label2: TLabel;
    Label3: TLabel;
    Label4: TLabel;
    Edithx: TEdit;
    Edithy: TEdit;
    Edithz: TEdit;
    Labelunion1: TLabel;
    Labelunion2: TLabel;
    Labelunion3: TLabel;
    Panelhxyz2: TPanel;
    Label5: TLabel;
    Label6: TLabel;
    Label7: TLabel;
    Label8: TLabel;
    Label9: TLabel;
    Label10: TLabel;
    Edit1: TEdit;
    Edit2: TEdit;
    Edit3: TEdit;
    Edit4: TEdit;
    Edit5: TEdit;
    Edit6: TEdit;
    Label11: TLabel;
    Label12: TLabel;
    Label13: TLabel;
    Label14: TLabel;
    Label15: TLabel;
    Label16: TLabel;
    procedure ButtonNextClick(Sender: TObject);
    procedure RadioGroup1Click(Sender: TObject);
  private
    { Private declarations }
  public
    { Public declarations }
  end;

var
  Formcabinetindent: TFormcabinetindent;

implementation

{$R *.dfm}

uses UnitVariables;

procedure TFormcabinetindent.ButtonNextClick(Sender: TObject);
var
   bOk : Boolean;
   r : Real;
   s : String;
   i : Integer;
begin
    case RadioGroup1.ItemIndex of
      0 : begin
             // без отступа.
              Close;
          end;
      1 : begin
             // с отступом.
             bOk:=true;
             s:=Trim(Edithonly.Text);
              if (FormatSettings.DecimalSeparator='.') then
              begin
                  for i:=1 to length(s) do
                  begin
                     if (s[i]=',') then s[i]:='.';
                  end;
                  Edithonly.Text:=Trim(s);
               end;
               if (FormatSettings.DecimalSeparator=',') then
               begin
                  for i:=1 to length(s) do
                  begin
                     if (s[i]='.') then s[i]:=',';
                  end;
                  Edithonly.Text:=Trim(s);
               end;
                if bOk then r:=FormVariables.my_real_convert(s,bOk);
                if (bOk) then
                begin
                   Edithonly.Color:=clwhite;
                   Close;
                end
                else
                begin
                   Edithonly.Color:=clred;
                end;
          end;
          2 : begin
              // отступ индивидуальный в каждом координатном направлении.
               bOk:=true;
              s:=Trim(Edithx.Text);
              if (FormatSettings.DecimalSeparator='.') then
              begin
                  for i:=1 to length(s) do
                  begin
                     if (s[i]=',') then s[i]:='.';
                  end;
                  Edithx.Text:=Trim(s);
               end;
               if (FormatSettings.DecimalSeparator=',') then
               begin
                  for i:=1 to length(s) do
                  begin
                     if (s[i]='.') then s[i]:=',';
                  end;
                  Edithx.Text:=Trim(s);
               end;
                if bOk then r:=FormVariables.my_real_convert(s,bOk);
                if (bOk) then
                begin
                   Edithx.Color:=clwhite;

                     bOk:=true;
                     s:=Trim(Edithy.Text);
                     if (FormatSettings.DecimalSeparator='.') then
                     begin
                        for i:=1 to length(s) do
                        begin
                            if (s[i]=',') then s[i]:='.';
                        end;
                        Edithy.Text:=Trim(s);
                     end;
                     if (FormatSettings.DecimalSeparator=',') then
                     begin
                        for i:=1 to length(s) do
                        begin
                           if (s[i]='.') then s[i]:=',';
                        end;
                        Edithy.Text:=Trim(s);
                     end;
                     if bOk then r:=FormVariables.my_real_convert(s,bOk);
                     if (bOk) then
                     begin
                        Edithy.Color:=clwhite;


                        bOk:=true;
                        s:=Trim(Edithz.Text);
                        if (FormatSettings.DecimalSeparator='.') then
                        begin
                           for i:=1 to length(s) do
                           begin
                               if (s[i]=',') then s[i]:='.';
                           end;
                           Edithz.Text:=Trim(s);
                        end;
                        if (FormatSettings.DecimalSeparator=',') then
                        begin
                           for i:=1 to length(s) do
                           begin
                              if (s[i]='.') then s[i]:=',';
                           end;
                           Edithz.Text:=Trim(s);
                        end;
                        if bOk then r:=FormVariables.my_real_convert(s,bOk);
                        if (bOk) then
                        begin
                           Edithz.Color:=clwhite;




                            Close;
                        end
                          else
                        begin
                           Edithz.Color:=clred;
                        end;




                     end
                      else
                     begin
                        Edithy.Color:=clred;
                     end;



                end
                else
                begin
                   Edithx.Color:=clred;
                end;
          end;
          3 : begin

              // отступ индивидуальный в каждом координатном направлении.
               bOk:=true;
              s:=Trim(Edit1.Text);
              if (FormatSettings.DecimalSeparator='.') then
              begin
                  for i:=1 to length(s) do
                  begin
                     if (s[i]=',') then s[i]:='.';
                  end;
                  Edit1.Text:=Trim(s);
               end;
               if (FormatSettings.DecimalSeparator=',') then
               begin
                  for i:=1 to length(s) do
                  begin
                     if (s[i]='.') then s[i]:=',';
                  end;
                  Edit1.Text:=Trim(s);
               end;
                if bOk then r:=FormVariables.my_real_convert(s,bOk);
                if (bOk) then
                begin
                   Edit1.Color:=clwhite;

                     bOk:=true;
                     s:=Trim(Edit2.Text);
                     if (FormatSettings.DecimalSeparator='.') then
                     begin
                        for i:=1 to length(s) do
                        begin
                            if (s[i]=',') then s[i]:='.';
                        end;
                        Edit2.Text:=Trim(s);
                     end;
                     if (FormatSettings.DecimalSeparator=',') then
                     begin
                        for i:=1 to length(s) do
                        begin
                           if (s[i]='.') then s[i]:=',';
                        end;
                        Edit2.Text:=Trim(s);
                     end;
                     if bOk then r:=FormVariables.my_real_convert(s,bOk);
                     if (bOk) then
                     begin
                        Edit2.Color:=clwhite;


                        bOk:=true;
                        s:=Trim(Edit3.Text);
                        if (FormatSettings.DecimalSeparator='.') then
                        begin
                           for i:=1 to length(s) do
                           begin
                               if (s[i]=',') then s[i]:='.';
                           end;
                           Edit3.Text:=Trim(s);
                        end;
                        if (FormatSettings.DecimalSeparator=',') then
                        begin
                           for i:=1 to length(s) do
                           begin
                              if (s[i]='.') then s[i]:=',';
                           end;
                           Edit3.Text:=Trim(s);
                        end;
                        if bOk then r:=FormVariables.my_real_convert(s,bOk);
                        if (bOk) then
                        begin
                           Edit3.Color:=clwhite;



                           // отступ индивидуальный в каждом координатном направлении.
                           bOk:=true;
                           s:=Trim(Edit4.Text);
                           if (FormatSettings.DecimalSeparator='.') then
                           begin
                              for i:=1 to length(s) do
                               begin
                                  if (s[i]=',') then s[i]:='.';
                               end;
                               Edit4.Text:=Trim(s);
                           end;
                           if (FormatSettings.DecimalSeparator=',') then
                           begin
                               for i:=1 to length(s) do
                               begin
                                  if (s[i]='.') then s[i]:=',';
                               end;
                               Edit4.Text:=Trim(s);
                           end;
                           if bOk then r:=FormVariables.my_real_convert(s,bOk);
                           if (bOk) then
                           begin
                              Edit4.Color:=clwhite;

                              bOk:=true;
                              s:=Trim(Edit5.Text);
                              if (FormatSettings.DecimalSeparator='.') then
                              begin
                                 for i:=1 to length(s) do
                                 begin
                                    if (s[i]=',') then s[i]:='.';
                                 end;
                                 Edit5.Text:=Trim(s);
                              end;
                              if (FormatSettings.DecimalSeparator=',') then
                              begin
                                 for i:=1 to length(s) do
                                 begin
                                    if (s[i]='.') then s[i]:=',';
                                 end;
                                 Edit5.Text:=Trim(s);
                              end;
                             if bOk then r:=FormVariables.my_real_convert(s,bOk);
                             if (bOk) then
                             begin
                                Edit5.Color:=clwhite;


                                bOk:=true;
                                s:=Trim(Edit6.Text);
                                if (FormatSettings.DecimalSeparator='.') then
                                begin
                                   for i:=1 to length(s) do
                                   begin
                                      if (s[i]=',') then s[i]:='.';
                                   end;
                                  Edit6.Text:=Trim(s);
                                end;
                                if (FormatSettings.DecimalSeparator=',') then
                                begin
                                   for i:=1 to length(s) do
                                   begin
                                      if (s[i]='.') then s[i]:=',';
                                   end;
                                   Edit6.Text:=Trim(s);
                                end;
                                if bOk then r:=FormVariables.my_real_convert(s,bOk);
                                if (bOk) then
                                begin
                                   Edit6.Color:=clwhite;




                                   Close;
                                end
                                 else
                                begin
                                   Edit6.Color:=clred;
                                end;




                             end
                              else
                             begin
                                Edit5.Color:=clred;
                             end;



                           end
                            else
                           begin
                              Edit4.Color:=clred;
                           end;




                            Close;
                        end
                          else
                        begin
                           Edit3.Color:=clred;
                        end;




                     end
                      else
                     begin
                        Edit2.Color:=clred;
                     end;



                end
                else
                begin
                   Edit1.Color:=clred;
                end;

          end;
   end;

end;

// смена типа отстпа.
procedure TFormcabinetindent.RadioGroup1Click(Sender: TObject);
begin
   case RadioGroup1.ItemIndex of
      0 : begin
             // без отступа.
             PanelH_only.Visible:=false;
             Panelhxyz.Visible:=false;
             Panelhxyz2.Visible:=false;
          end;
      1 : begin
             // с отступом.
             PanelH_only.Visible:=true;
             Panelhxyz.Visible:=false;
             Panelhxyz2.Visible:=false;
          end;
      2 : begin
             // hx;hy;hz;
             PanelH_only.Visible:=false;
             Panelhxyz.Visible:=true;
             Panelhxyz2.Visible:=false;
          end;
      3 : begin
             // hxmin; hxmax; hymin; hymax; hzmin; hzmax;
             PanelH_only.Visible:=false;
             Panelhxyz.Visible:=false;
             Panelhxyz2.Visible:=true;
          end;
   end;
end;

end.
