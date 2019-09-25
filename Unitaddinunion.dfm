object Form_formirateunion: TForm_formirateunion
  Left = 291
  Top = 126
  AutoSize = True
  Caption = 'Add in union'
  ClientHeight = 89
  ClientWidth = 194
  Color = clMoneyGreen
  Font.Charset = DEFAULT_CHARSET
  Font.Color = clWindowText
  Font.Height = -11
  Font.Name = 'Tahoma'
  Font.Style = []
  OldCreateOrder = False
  OnCreate = FormCreate
  PixelsPerInch = 96
  TextHeight = 13
  object grp1: TGroupBox
    Left = 0
    Top = 0
    Width = 185
    Height = 89
    Caption = 'Please, select union:'
    Color = clMoneyGreen
    ParentColor = False
    TabOrder = 0
    object cbbselectnameunion: TComboBox
      Left = 8
      Top = 24
      Width = 169
      Height = 21
      TabOrder = 0
    end
    object btnApply: TButton
      Left = 96
      Top = 56
      Width = 75
      Height = 25
      Caption = 'Add'
      TabOrder = 1
      OnClick = btnApplyClick
    end
  end
  object ApplicationEvents1: TApplicationEvents
    OnMessage = ApplicationEvents1Message
    Left = 8
    Top = 56
  end
end
