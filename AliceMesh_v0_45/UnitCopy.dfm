object FormCopyObject: TFormCopyObject
  Left = 336
  Top = 183
  AutoSize = True
  Caption = 'Copy object'
  ClientHeight = 195
  ClientWidth = 251
  Color = clMoneyGreen
  CustomTitleBar.CaptionAlignment = taCenter
  Font.Charset = DEFAULT_CHARSET
  Font.Color = clWindowText
  Font.Height = -11
  Font.Name = 'MS Sans Serif'
  Font.Style = []
  OldCreateOrder = False
  OnClose = FormClose
  OnCreate = FormCreate
  PixelsPerInch = 96
  TextHeight = 13
  object PanelMain: TPanel
    Left = 0
    Top = 0
    Width = 251
    Height = 195
    AutoSize = True
    Color = clMoneyGreen
    TabOrder = 0
    object Labelnc: TLabel
      Left = 1
      Top = 9
      Width = 83
      Height = 13
      Caption = 'Number of copies'
    end
    object Editnc: TEdit
      Left = 97
      Top = 1
      Width = 73
      Height = 21
      TabOrder = 0
      Text = '1'
    end
    object GroupBoxTranslate: TGroupBox
      Left = 1
      Top = 28
      Width = 249
      Height = 121
      Caption = 'Translate'
      Color = clMoneyGreen
      ParentBackground = False
      ParentColor = False
      TabOrder = 1
      object Labelx: TLabel
        Left = 8
        Top = 24
        Width = 33
        Height = 13
        Caption = 'Xoffset'
      end
      object Labely: TLabel
        Left = 8
        Top = 56
        Width = 33
        Height = 13
        Caption = 'Yoffset'
      end
      object Labelz: TLabel
        Left = 8
        Top = 88
        Width = 33
        Height = 13
        Caption = 'Zoffset'
      end
      object Label1: TLabel
        Left = 176
        Top = 24
        Width = 32
        Height = 13
        Caption = 'Label1'
      end
      object Label2: TLabel
        Left = 176
        Top = 56
        Width = 32
        Height = 13
        Caption = 'Label2'
      end
      object Label3: TLabel
        Left = 176
        Top = 88
        Width = 32
        Height = 13
        Caption = 'Label3'
      end
      object EditX: TEdit
        Left = 56
        Top = 16
        Width = 113
        Height = 21
        TabOrder = 0
        Text = '0.0'
      end
      object EditY: TEdit
        Left = 56
        Top = 48
        Width = 113
        Height = 21
        TabOrder = 1
        Text = '0.0'
      end
      object EditZ: TEdit
        Left = 56
        Top = 80
        Width = 113
        Height = 21
        TabOrder = 2
        Text = '0.0'
      end
    end
    object ButtonApply: TButton
      Left = 97
      Top = 169
      Width = 75
      Height = 25
      Caption = 'Apply'
      TabOrder = 2
      OnClick = ButtonApplyClick
    end
    object CheckBoxRotate: TCheckBox
      Left = 1
      Top = 173
      Width = 64
      Height = 17
      Caption = 'Rotate'
      TabOrder = 3
      Visible = False
    end
  end
end
